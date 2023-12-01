import pytest
import numpy as np
import tensorflow as tf

from pytest import approx

from hmsc.updaters.updateZ import updateZ


def input_values(*, distr_case, x_ndim, seed, dtype):
    assert distr_case in (0, 1, 2, 3)
    assert x_ndim in (2, 3)
    # Add test Y including nans

    rng = np.random.default_rng(seed=seed)
    ns = 7
    ny = 5
    nb = 11
    nr = 2
    ne = 3

    distr = np.array(
        [[1, 1],
         [2, 1],
         [3, 1],
         [1, 1],
         [3, 1],
         [2, 1],
         [1, 1]]
    )

    if distr_case in (1, 2, 3):
        distr[:, 0] = distr_case

    Pi = np.array(
        [[0, 2],
         [1, 2],
         [2, 2],
         [0, 0],
         [1, 0]]
    )

    ns3 = np.sum(distr[:, 0] == 3)

    params = dict(
        Z=tf.convert_to_tensor(rng.random(size=(ny, ns)), dtype=dtype),
        Beta=tf.convert_to_tensor(rng.random(size=(nb, ns)), dtype=dtype),
        sigma=tf.convert_to_tensor(rng.random(size=ns), dtype=dtype),
        Xeff=tf.convert_to_tensor(rng.random(size=(ny, nb)), dtype=dtype),
        poisson_omega=tf.convert_to_tensor(rng.random(size=(ny, ns3)), dtype=dtype),
    )
    data = dict(
        Y=tf.convert_to_tensor(rng.integers(2, size=(ny, ns)), dtype=dtype),
        Pi=tf.convert_to_tensor(Pi, dtype=tf.int32),
        distr=tf.convert_to_tensor(distr, dtype=tf.int32),
    )

    assert data['distr'].shape == (ns, 2)
    assert data['Pi'].shape == (ny, nr)

    EtaList = [
        tf.convert_to_tensor(rng.random(ne * 2).reshape(ne, 2), dtype=dtype),
        tf.convert_to_tensor(rng.random(ne * 2).reshape(ne, 2), dtype=dtype),
    ]
    LambdaList = [
        tf.convert_to_tensor(rng.random(2 * ns).reshape(2, ns), dtype=dtype),
        tf.convert_to_tensor(rng.random(2 * ns).reshape(2, ns), dtype=dtype),
    ]
    rLHyperparams = [
        dict(xDim=0),
        dict(xDim=0),
    ]
    assert len(EtaList) == len(LambdaList) == len(rLHyperparams) == nr

    if x_ndim == 3:
        params['Xeff'] = tf.convert_to_tensor(rng.random(size=(ns, ny, nb)), dtype=dtype)

    params["Eta"] = EtaList
    params["Lambda"] = LambdaList

    return params, data, rLHyperparams


def reference_values(*,
        distr_case,
        x_ndim,
        poisson_preupdate_z,
        poisson_marginalize_z,
        truncated_normal_library,
        seed,
    ):
    tnlib = truncated_normal_library
    if seed != 42:
        raise ValueError('reference missing')

    if x_ndim == 3 and distr_case == 0 and tnlib == 'tf':
        Z = \
[[ 0.        ,  3.53200434,  1.75342901,  1.        ,  0.13864651, -0.0499102 ,  1.        ],
 [ 0.        ,  2.64775599,  1.43855754,  0.        ,  0.87347456, -0.20061524,  1.        ],
 [ 1.        , -0.01429325,  1.72133414,  0.        ,  0.64247312,  3.83834427,  1.        ],
 [ 0.        ,  2.16941234,  1.57816434,  0.        ,  0.14337646,  2.3323683 ,  1.        ],
 [ 0.        , -0.06552546,  1.63378067,  0.        ,  0.29615809,  2.05435849,  0.        ]]
    elif distr_case == 0 and tnlib == 'tf':
        Z = \
[[ 0.        ,  2.53617218,  1.40754338,  1.        ,  0.13129104, -0.07285226,  1.        ],
 [ 0.        ,  2.97700928,  1.35676311,  0.        ,  0.89406079, -0.16933972,  1.        ],
 [ 1.        , -0.01509904,  1.77882987,  0.        ,  0.6501287 ,  3.88307341,  1.        ],
 [ 0.        ,  2.47942512,  1.68541477,  0.        ,  0.11587453,  3.06711182,  1.        ],
 [ 0.        , -0.06391573,  1.44970694,  0.        ,  0.32676159,  2.30020615,  0.        ]]
    elif distr_case == 0 and tnlib == 'tfd':
        Z = \
[[ 0.        ,  2.17921306,  1.40754338,  1.        ,  0.13129104, -0.05694204,  1.        ],
 [ 0.        ,  2.7920181 ,  1.35676311,  0.        ,  0.89406079, -0.13361183,  1.        ],
 [ 1.        , -0.0166318 ,  1.77882987,  0.        ,  0.6501287 ,  4.70081388,  1.        ],
 [ 0.        ,  4.68685993,  1.68541477,  0.        ,  0.11587453,  4.4379023 ,  1.        ],
 [ 0.        , -0.1906678 ,  1.44970694,  0.        ,  0.32676159,  2.69563436,  0.        ]]
    elif distr_case == 1:
        Z = \
[[0., 0., 1., 0., 0., 0., 1.],
 [1., 0., 1., 0., 1., 1., 1.],
 [0., 1., 0., 1., 0., 0., 0.],
 [1., 0., 1., 0., 0., 1., 0.],
 [1., 1., 0., 0., 0., 1., 1.]]
    elif distr_case == 2 and tnlib == 'tf':
        Z = \
[[-2.16731687e-02, -7.26993294e-02,  2.97247856e+00, -2.56566343e-01, -2.82493186e-02, -2.15090045e-01,  1.82051501e+00],
 [ 3.44111121e+00, -8.08565606e-02,  3.02281863e+00, -1.12269164e-01,  3.85723288e+00,  3.10194142e+00,  1.86781389e+00],
 [-1.97487684e-02,  3.19087864e+00, -7.81100111e-03,  3.84275918e+00, -1.37848714e-01, -1.85120295e-02, -9.35641902e-02],
 [ 2.98396344e+00, -2.44607084e-03,  4.53406253e+00, -2.99807014e-01, -1.35385546e-01,  2.85929925e+00, -1.04299383e-01],
 [ 3.86621611e+00,  4.19202896e+00, -1.22456847e-02, -1.25848590e-01, -1.53532906e-01,  2.53178412e+00,  2.94148879e+00]]
    elif distr_case == 3:
        Z = \
[[1.0157138 , 0.4100475 , 1.3975398 , 0.45718054, 0.14703785, 0.89889734, 0.77702612],
 [0.93864013, 0.16124331, 1.27282378, 0.28755889, 0.83547283, 0.71256083, 0.98335849],
 [0.93121183, 0.47809308, 1.93095805, 0.29002695, 1.18488449, 0.65230171, 0.92880115],
 [0.42809156, 1.1224092 , 1.78635417, 0.70666318, 0.18722986, 0.63402328, 0.11470677],
 [0.5574293 , 0.85028616, 1.43130618, 0.87880138, 0.44498424, 0.25420793, 0.84656094]]

    iD = \
[[10.68268147,  2.98066761, 32.00141647,  1.36279137,  1.73806681,  1.93189325,  5.35606031],
 [10.68268147,  2.98066761, 32.00141647,  1.36279137,  1.73806681,  1.93189325,  5.35606031],
 [10.68268147,  2.98066761, 32.00141647,  1.36279137,  1.73806681,  1.93189325,  5.35606031],
 [10.68268147,  2.98066761, 32.00141647,  1.36279137,  1.73806681,  1.93189325,  5.35606031],
 [10.68268147,  2.98066761, 32.00141647,  1.36279137,  1.73806681,  1.93189325,  5.35606031]]


    if distr_case == 0:
        omega = \
[[82.69564925, 73.52446685],
 [75.73119083, 80.31543575],
 [77.77786621, 79.68186375],
 [81.55596303, 73.14161346],
 [80.27240777, 75.35136998]]
    elif distr_case in [1, 2]:
        omega = \
[]
    elif distr_case in [3]:
        omega = \
[[81.70836675, 75.63906919, 81.970256  , 79.53973508, 73.76147381, 84.78603609, 80.57772435],
 [81.23612849, 73.43414908, 76.19355286, 74.56216626, 80.7187692 , 78.19513098, 83.49581729],
 [77.96440037, 75.55425842, 79.08140616, 74.00632178, 84.13732389, 78.59819836, 83.34528905],
 [75.25213129, 84.97194624, 82.77510123, 80.53394158, 73.69114804, 77.7371134 , 72.00083597],
 [72.77825139, 82.103995  , 79.33289536, 82.3816733 , 76.11819717, 75.78188646, 81.89521818]]


    if distr_case == 3 and poisson_preupdate_z and poisson_marginalize_z:
        Z = \
[[ -43.42004878, -162.94508265,  -12.30981522, -207.30316425, -185.01985086, -181.6067665 ,  -76.38012609],
 [ -39.70457435, -139.31378297,  -12.01302183, -230.3347691 , -204.78232546, -214.0327703 ,  -88.12797169],
 [ -39.94141362, -149.98152148,  -10.09754038, -306.65202511, -262.82041279, -175.75076706,  -85.74530644],
 [ -43.48429453, -134.21369115,  -10.70082485, -317.11661209, -240.39275032, -179.1056778 ,  -74.84870778],
 [ -41.22543472, -154.71900044,  -12.13989044, -310.89598669, -262.18071751, -187.19134434,  -86.74641565]]
        iD = \
[[ 5.14760595,  1.48029499, 14.34262759,  0.86011242,  1.04211029,  1.11727675,  2.83061468],
 [ 5.35234156,  1.59177087, 14.46590774,  0.82762748,  1.00084487,  1.0422008 ,  2.65275836],
 [ 5.33613414,  1.54018213, 15.32350804,  0.73480375,  0.89658099,  1.1321058 ,  2.68809604],
 [ 5.14420346,  1.61806871, 15.04533678,  0.72367448,  0.93462059,  1.12356242,  2.85403605],
 [ 5.26659959,  1.51728399, 14.41295148,  0.73024921,  0.89761165,  1.10349244,  2.67370574]]
        omega = \
[[ 9.93486621,  2.94078101, 25.99183908,  2.3318139 ,  2.60254379,  2.64966325,  6.00327433],
 [10.7267756 ,  3.41604941, 26.39954993,  2.10754797,  2.35958151,  2.26305246,  5.25591813],
 [10.66187537,  3.18696096, 29.40260551,  1.59459226,  1.85186443,  2.73461098,  5.39647574],
 [ 9.92220007,  3.53950454, 28.39524802,  1.54309382,  2.02183169,  2.68528994,  6.10960872],
 [10.38784257,  3.09045367, 26.22371328,  1.57329803,  1.85626681,  2.5734277 ,  5.33879052]]

    if distr_case == 3 and not poisson_preupdate_z and poisson_marginalize_z:
        Z = \
[[ 0.78843091,  0.30402515,  0.8140815 ,  0.62787517,  0.13592746,  1.01645498,  0.70256639],
 [ 0.75285832,  0.10574278,  0.35208251,  0.20194153,  0.71960344,  0.51349526,  0.92541968],
 [ 0.50098504,  0.28999461,  0.58515653,  0.1515757 ,  0.97103211,  0.55264778,  0.90861597],
 [ 0.26342519,  1.02934455,  0.86729151,  0.6991929 ,  0.12267982,  0.48225299, -0.02966417],
 [ 0.03757098,  0.82400747,  0.61150203,  0.83844423,  0.3455925 ,  0.31647043,  0.80239271]]
        iD = \
[[ 9.4475003 ,  2.86766317, 23.01593232,  1.33983535,  1.69805496,  1.88885473,  5.02222907],
 [ 9.44115448,  2.86440247, 22.5361829 ,  1.33833037,  1.70143096,  1.88531456,  5.03319283],
 [ 9.39533302,  2.86754113, 22.7822534 ,  1.33814998,  1.70288937,  1.8855477 ,  5.03264491],
 [ 9.35470181,  2.87965436, 23.07894107,  1.34011403,  1.69801765,  1.88504678,  4.98521578],
 [ 9.31533894,  2.8762495 , 22.80307825,  1.34061438,  1.69926612,  1.88386815,  5.02726992]]
        omega = \
[[81.70836675, 75.63906919, 81.970256  , 79.53973508, 73.76147381, 84.78603609, 80.57772435],
 [81.23612849, 73.43414908, 76.19355286, 74.56216626, 80.7187692 , 78.19513098, 83.49581729],
 [77.96440037, 75.55425842, 79.08140616, 74.00632178, 84.13732389, 78.59819836, 83.34528905],
 [75.25213129, 84.97194624, 82.77510123, 80.53394158, 73.69114804, 77.7371134 , 72.00083597],
 [72.77825139, 82.103995  , 79.33289536, 82.3816733 , 76.11819717, 75.78188646, 81.89521818]]

    if distr_case == 3 and poisson_preupdate_z and not poisson_marginalize_z:
        Z = \
[[ -19.66038667,  -79.88297817,   -3.67414683, -128.95783388, -109.20056722, -104.02445539,  -38.84796085],
 [ -18.32891211,  -73.44637301,   -3.49158543, -137.98937781, -116.23370342, -114.52806157,  -42.07514868],
 [ -18.16217468,  -75.17924591,   -2.60915876, -163.37094787, -134.36329774, -100.82970349,  -41.21746715],
 [ -19.51294963,  -70.89813222,   -2.34664345, -167.55925337, -128.03929231, -102.58689855,  -37.8535931 ],
 [ -18.74517299,  -76.75401524,   -3.52758133, -164.66722095, -133.24234651, -105.24438911,  -41.42424892]]
        iD = \
[[10.68268147,  2.98066761, 32.00141647,  1.36279137,  1.73806681,  1.93189325,  5.35606031],
 [10.68268147,  2.98066761, 32.00141647,  1.36279137,  1.73806681,  1.93189325,  5.35606031],
 [10.68268147,  2.98066761, 32.00141647,  1.36279137,  1.73806681,  1.93189325,  5.35606031],
 [10.68268147,  2.98066761, 32.00141647,  1.36279137,  1.73806681,  1.93189325,  5.35606031],
 [10.68268147,  2.98066761, 32.00141647,  1.36279137,  1.73806681,  1.93189325,  5.35606031]]
        omega = \
[[ 9.93486621,  2.94078101, 25.99183908,  2.3318139 ,  2.60254379,  2.64966325,  6.00327433],
 [10.7267756 ,  3.41604941, 26.39954993,  2.10754797,  2.35958151,  2.26305246,  5.25591813],
 [10.66187537,  3.18696096, 29.40260551,  1.59459226,  1.85186443,  2.73461098,  5.39647574],
 [ 9.92220007,  3.53950454, 28.39524802,  1.54309382,  2.02183169,  2.68528994,  6.10960872],
 [10.38784257,  3.09045367, 26.22371328,  1.57329803,  1.85626681,  2.5734277 ,  5.33879052]]


    for var in ['Z', 'iD', 'omega']:
        if var not in locals():
            raise ValueError(f'reference for {var} missing')

    return np.array(Z), np.array(iD), np.array(omega)


def run_test(distr_case=0, x_ndim=2, tnlib='tf',
             poisson_preupdate_z=False,
             poisson_marginalize_z=False,
             seed=42, dtype=np.float64):
    input_kwargs = dict(
        distr_case=distr_case,
        x_ndim=x_ndim,
        seed=seed,
    )
    func_kwargs = dict(
        poisson_preupdate_z=poisson_preupdate_z,
        poisson_marginalize_z=poisson_marginalize_z,
        truncated_normal_library=tnlib,
        seed=seed,
    )
    ref_kwargs = input_kwargs.copy()
    ref_kwargs.update(func_kwargs)

    # Prepare input matrices
    params, data, rLHyperparams = input_values(**input_kwargs, dtype=dtype)

    # Calculate
    Z, iD, omega = updateZ(
        params, data, rLHyperparams,
        dtype=dtype, **func_kwargs,
        )
    Z = Z.numpy()
    iD = iD.numpy()
    omega = omega.numpy()

    assert Z.shape == params['Z'].shape
    assert iD.shape == params['Z'].shape
    assert omega.shape == params['poisson_omega'].shape

    # Print values
    for name, array in (('Z', Z),
                        ('iD', iD),
                        ('omega', omega)):
        print(f'{name} = \\')
        print(np.array2string(array, separator=', ', max_line_width=200))

    # Test against reference
    ref_Z, ref_iD, ref_omega = reference_values(**ref_kwargs)
    assert Z == approx(ref_Z)
    assert iD == approx(ref_iD)
    if omega.size > 0 or ref_omega.size > 0:
        assert omega == approx(ref_omega)


def test_simple():
    run_test()


@pytest.mark.parametrize("tnlib", ['tf', 'tfd'])
def test_tnlib(tnlib):
    run_test(tnlib=tnlib)


@pytest.mark.parametrize("distr_case", [0, 1, 2, 3])
def test_distr_case(distr_case):
    run_test(distr_case=distr_case)


@pytest.mark.parametrize("x_ndim", [2, 3])
def test_x_ndim(x_ndim):
    run_test(x_ndim=x_ndim)


@pytest.mark.parametrize("preupdate", [True, False])
@pytest.mark.parametrize("marginalize", [True, False])
def test_poisson_flags(preupdate, marginalize):
    run_test(
        distr_case=3,
        poisson_preupdate_z=preupdate,
        poisson_marginalize_z=marginalize,
    )
