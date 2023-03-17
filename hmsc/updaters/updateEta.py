import numpy as np
import tensorflow as tf

tfla, tfm, tfr, tfs = tf.linalg, tf.math, tf.random, tf.sparse

#from hmsc.utils.tflautils import kron, scipy_cholesky, tf_sparse_matmul, tf_sparse_cholesky, scipy_sparse_solve_triangular, convert_sparse_tensor_to_sparse_csc_matrix
from hmsc.utils.tflautils import tf_sparse_matmul, tf_sparse_cholesky, scipy_sparse_solve_triangular

from scipy.sparse.linalg import splu, spsolve_triangular
from scipy.sparse import csc_matrix, coo_matrix, block_diag, kron

def updateEta(params, data, modelDims, rLHyperparams, dtype=np.float64):
    """Update conditional updater(s):
    Z - site loadings.

    Parameters
    ----------
    params : dict
        The initial value of the model parameter(s):
        Z - latent variables
        Beta - species niches
        Eta - site loadings
        Lambda - species loadings
        Alpha - scale of site loadings (eta's prior)
        sigma - residual variance
        X - environmental data
        Pi - study design
        iWg - ??
        sDim - ??
    """
    Z = params["Z"]
    Beta = params["Beta"]
    sigma = params["sigma"]
    LambdaList = params["Lambda"]
    EtaList = params["Eta"]
    AlphaIndList = params["AlphaInd"]
    Pi = data["Pi"]
    X = data["X"]
    Y = data["Y"]
    ny = modelDims["ny"]
    ns = modelDims["ns"]
    nr = modelDims["nr"]
    npVec = modelDims["np"]
    
    iD = tf.ones_like(Z) * sigma**-2
    
    if isinstance(X, list):
        LFix = tf.einsum("jik,kj->ij", tf.stack(X), Beta)
    else:
        LFix = tf.matmul(X, Beta)

    LRanLevelList = [None] * nr
    for r, (Eta, Lambda) in enumerate(zip(EtaList, LambdaList)):
        LRanLevelList[r] = tf.matmul(tf.gather(Eta, Pi[:,r]), Lambda)

    EtaListNew = [None] * nr
    for r, (Eta, Lambda, AlphaInd, rLPar) in enumerate(zip(EtaList, LambdaList, AlphaIndList, rLHyperparams)):
        nf = tf.cast(tf.shape(Lambda)[-2], tf.int64)
        if nf > 0:
            S = Z - LFix - sum([LRanLevelList[rInd] for rInd in np.setdiff1d(np.arange(nr), r)])
            mu0 = tf.scatter_nd(Pi[:,r,None], tf.matmul(iD * S, Lambda, transpose_b=True), [npVec[r],nf])
            LamInvSigLam = tf.scatter_nd(Pi[:,r,None], tf.einsum("hj,ij,kj->ihk", Lambda, iD, Lambda), [npVec[r],nf,nf])
            
            if rLPar["sDim"] == 0:
                EtaListNew[r] = modelNonSpatial(LamInvSigLam, mu0, npVec[r], nf, dtype)
            else:
                if rLPar["spatialMethod"] == "Full":
                    EtaListNew[r] = modelSpatialFull(LamInvSigLam, mu0, AlphaInd, rLPar["iWg"], npVec[r], nf)
                elif rLPar["spatialMethod"] == "GPP":
                    EtaListNew[r] = modelSpatialGPP(LamInvSigLam, mu0, AlphaInd, rLPar["Fg"], rLPar["idDg"], rLPar["idDW12g"], rLPar["nK"], npVec[r], nf)
                    # EtaListNew[r] = modelSpatialGPP(Lambda, AlphaInd, iSigma, mu0, rLPar["Fg"], rLPar["idDg"], rLPar["idDW12g"], npVec[r], nf, rLPar["nK"])
                elif rLPar["spatialMethod"] == "NNGP":                
                    modelSpatialNNGP_local = lambda LamInvSigLam, mu0, Alpha, nf: modelSpatialNNGP_scipy(LamInvSigLam, mu0, Alpha, rLPar["iWList_csc"], npVec[r], nf)
                    # EtaListNew[r] = modelSpatialNNGP_local(LamInvSigLam, mu0, AlphaInd, nf)
                    Eta = tf.numpy_function(modelSpatialNNGP_local, [LamInvSigLam, mu0, AlphaInd, nf], dtype)
                    EtaListNew[r] = tf.ensure_shape(Eta, [npVec[r], None])              
            
            LRanLevelList[r] = tf.matmul(tf.gather(EtaListNew[r], Pi[:,r]), Lambda)
        else:
            EtaListNew[r] = Eta

    return EtaListNew


def modelNonSpatial(LamInvSigLam, mu0, np, nf, dtype=np.float64):
    iV = tf.eye(nf, dtype=dtype) + LamInvSigLam
    LiV = tfla.cholesky(iV)
    mu1 = tfla.triangular_solve(LiV, tf.expand_dims(mu0, -1))
    Eta = tf.squeeze(tfla.triangular_solve(LiV, mu1 + tfr.normal([np,nf,1], dtype=dtype), adjoint=True), -1)
    return Eta


def modelSpatialFull(LamInvSigLam, mu0, AlphaInd, iWg, np, nf, dtype=np.float64): 
    #TODO a lot of unnecessary tanspositions - rework if considerably affects perfomance
    iWs = tf.reshape(tf.transpose(tfla.diag(tf.transpose(tf.gather(iWg, AlphaInd), [1,2,0])), [2,0,3,1]), [nf*np,nf*np])
    iUEta = iWs + tf.reshape(tf.transpose(tfla.diag(tf.transpose(LamInvSigLam, [1,2,0])), [0,2,1,3]), [nf*np,nf*np])
    LiUEta = tfla.cholesky(iUEta)
    mu1 = tfla.triangular_solve(LiUEta, tf.reshape(tf.transpose(mu0), [nf*np,1]))
    eta = tfla.triangular_solve(LiUEta, mu1 + tfr.normal([nf*np,1], dtype=dtype), adjoint=True)
    Eta = tf.transpose(tf.reshape(eta, [nf,np]))
    return Eta


def modelSpatialNNGP_scipy(LamInvSigLam, mu0, Alpha, iWList, nu, nf, dtype=np.float64):
    LamInvSigLam_bdiag = block_diag([LamInvSigLam[i] for i in range(nu)], dtype=dtype)
    dataList, colList, rowList = [None]*int(nf), [None]*int(nf), [None]*int(nf)
    for h, a in enumerate(Alpha):
      iW = coo_matrix(iWList[a])
      dataList[h] = iW.data
      colList[h] = iW.col + h*nu
      rowList[h] = iW.row + h*nu
    dataArray = np.concatenate(dataList)
    colArray = np.concatenate(colList)
    rowArray = np.concatenate(rowList)
    iUEta = csc_matrix(coo_matrix((dataArray,(colArray,rowArray)), [nu*nf,nu*nf])) + LamInvSigLam_bdiag
    # iWs = [kron(iWList[a], csc_matrix(coo_matrix(([1],([h],[h])), [nf,nf]))) for h, a in enumerate(Alpha)] #replaced with indices
    # iUEta = sum(iWs) + LamInvSigLam_bdiag
    LU_factor = splu(iUEta, "NATURAL", diag_pivot_thresh=0)
    L, U = LU_factor.L, LU_factor.U
    LiUEta = L.multiply(np.sqrt(U.diagonal()))
    mu1 = spsolve_triangular(LiUEta, np.reshape(mu0, [nu*nf]))
    eta = spsolve_triangular(LiUEta.transpose(), mu1 + np.random.normal(dtype(0), dtype(1), size=[nf*nu]), lower=False)
    Eta = np.reshape(eta, [nu,nf])
    return Eta


def modelSpatialGPP(LamInvSigLam, mu0, AlphaInd, Fg, idDg, idDW12g, nu, nf, nK, dtype=tf.float64):
    idD = tf.gather(idDg, AlphaInd, axis=1)
    Fmat = tf.reshape(tf.transpose(tfla.diag(tf.transpose(tf.gather(Fg, AlphaInd), [1,2,0])), [2,0,3,1]), [nf*nK,nf*nK])
    idD1W12 = tf.reshape(tf.transpose(tfla.diag(tf.transpose(tf.gather(idDW12g, AlphaInd), [1,2,0])), [2,0,3,1]), [nf*nK,nf*nu])


    B = tf.tile(LamSigLamT[None,:], [nu,1,1]) + tfla.diag(idD, num_rows = tf.cast(nf, tf.int32), num_cols = tf.cast(nf, tf.int32))
    LB = tfla.cholesky(B)
    iLB = tfla.inv(LB)

    LA = tf.reshape(tf.transpose(tfla.diag(tf.transpose(LB, [1,2,0])), [2,0,3,1]), [nf*nu,nf*nu])
    iLA = tf.reshape(tf.transpose(tfla.diag(tf.transpose(iLB, [1,2,0])), [2,0,3,1]), [nf*nu,nf*nu])

    iAidD1W12 = tf.einsum("ij,kj->ik", iLA, idD1W12)
    H = Fmat - tf.matmul(idD1W12, iAidD1W12)
    LH = tfla.cholesky(H)
    iLH = tfla.triangular_solve(LH, tf.eye(nK*nf, dtype=dtype), adjoint=True)

    mu1 = tf.matmul(iLA, tf.reshape(mu0,[-1])[:,None])
    tmp1 = tf.matmul(iAidD1W12, iLH)
    mu2 = tf.einsum("ij,kj,kl->il", tmp1, tmp1, tf.reshape(mu0,[-1])[:,None])

    etaR = tf.matmul(iLA, tfr.normal([nu*nf,1], dtype=dtype)) + tf.matmul(tmp1, tfr.normal([nK*nf,1], dtype=dtype))
    Eta = tf.reshape(mu1 + mu2 + etaR, shape=[nu,nf])

    return Eta


# Anis's NNGP version, not used currently, but contains potentially useful parts
def modelSpatialNNGP(S, iD, Pi, Lambda, AlphaInd, iWg, iSigma, ny, ns, nu, nf, dtype=np.float64):
    LamInvSigLam = tfs.from_dense(tf.scatter_nd(Pi, tf.einsum("hj,ij,kj->ihk", Lambda, iD, Lambda), [nu,nf,nf]))
    mu0 = tf.scatter_nd(Pi, tf.matmul(iD * S, Lambda, transpose_b=True), [nu,nf])

    mask_values = tf.constant([True, False])
    def tf_sparse_gather(A, Ind):
        def single_iter(x):
            mask = tf.where(tf.equal(A.indices[:, 0], tf.cast(Ind[x], dtype=tf.int64)), mask_values[0], mask_values[1])
            return tfs.to_dense(tfs.reduce_sum(tf.sparse.retain(A, mask), axis=0, keepdims=True, output_is_sparse=True))[-1, :, :]
        return tfs.from_dense(tf.map_fn(single_iter, tf.range(nf), fn_output_signature=dtype))

    def tf_sparse_block_diag(A, n, nx, ny):
        def tf_sparse_update_indices():
            x_stride = A.indices[:,0]*nx
            y_stride = A.indices[:,0]*ny
            return tf.stack([A.indices[:,1]+x_stride, A.indices[:,2]+y_stride], axis=1)
        return tfs.SparseTensor(indices=tf_sparse_update_indices(), values=A.values, dense_shape=(n*nx,n*ny))

    LamInvSigLam_bdiag = tf_sparse_block_diag(LamInvSigLam, nu, nf, nf)
    iWg_bdiag = tf_sparse_block_diag(tf_sparse_gather(iWg, AlphaInd), nf, nu, nu)
    iUEta = tfs.add(tfs.add(iWg_bdiag, LamInvSigLam_bdiag), tfs.eye(nu*nf, dtype=dtype))   

    LiUEta = tf_sparse_cholesky(iUEta)
    mu1 = tf.numpy_function(scipy_sparse_solve_triangular, [LiUEta.values, LiUEta.indices, [ny*ns, ny*ns], tf.reshape(tf.transpose(mu0), [nf * nu, 1])], dtype)
    Eta = tf.numpy_function(scipy_sparse_solve_triangular, [LiUEta.values, LiUEta.indices, [ny*ns, ny*ns], mu1 + tfr.normal([nf * nu, 1], dtype=dtype)], dtype)

    return tf.transpose(tf.reshape(Eta, [nf, nu]))