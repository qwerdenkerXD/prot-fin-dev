from scipy import signal
import numpy as np
from tools import *
from os import environ as env
import pandas as pd
from operator import gt, lt

WINDOW_TYPE = env.get("WINDOW_TYPE", ["boxcar", "triang", "blackman", "hamming", "hann", "bartlett", "flattop", "parzen", "bohman", "blackmanharris", "nuttall", "barthann", "cosine", "exponential", "tukey", "taylor", "lanczos"]\
                                     [0])

AMPL_QUANTILES = {
    10: {
        5: [ 0.844, 0.54312, 0.0482, 0.1085, 0.04832, 0.025],
        95: [ 1.576, 1.0551, 0.36524, 0.58613, 0.35894, 0.59],
        0.1: [ 0.607, 0.3849, 0.00222, 0.01523, 0.00222, 0.0],
        99.9: [ 1.908, 1.26526, 0.48708, 0.75536, 0.48234, 0.869],
        0.01: [ 0.525, 0.33408, 0.0, 0.00457, 0.0, 0.0],
        99.99: [ 2.03, 1.32497, 0.52196, 0.80759, 0.51932, 0.964],
        0.001: [ 0.469, 0.29125, 0.0, 4e-05, 0.0, 0.0],
        99.999: [ 2.1, 1.35915, 0.5426, 0.83649, 0.54057, 1.03],
        "std_dev": [ 0.22193044689998173, 0.15550996957329452, 0.09710811853360825, 0.1440028540709493, 0.09512672500832252, 0.1763286767067828],
        "mean": [ 1.201530310065123, 0.7917699436132187, 0.19102350789394956, 0.3386462044099786, 0.18884220292739026, 0.26475480031296555],
    },
    20: {
        5: [ 0.9265, 0.58593, 0.03295, 0.12298, 0.0336, 0.0614, 0.03466, 0.05002, 0.0334, 0.04578, 0.009],
        95: [ 1.47, 0.96363, 0.25633, 0.45731, 0.25985, 0.37825, 0.2688, 0.33932, 0.25623, 0.32337, 0.2875],
        0.1: [ 0.708, 0.44335, 0.0034, 0.02318, 0.00383, 0.00863, 0.00357, 0.00705, 0.0038, 0.00641, 0.0],
        99.9: [ 1.7525, 1.14938, 0.37494, 0.60134, 0.37683, 0.52166, 0.38326, 0.47933, 0.37019, 0.45956, 0.475],
        0.01: [ 0.6305, 0.39325, 0.0, 0.00746, 0.0, 0.00255, 0.0, 0.00223, 0.0, 0.00203, 0.0],
        99.99: [ 1.913, 1.23955, 0.42063, 0.6568, 0.42083, 0.57858, 0.4248, 0.53169, 0.41367, 0.51087, 0.5535],
        0.001: [ 0.5645, 0.35104, 0.0, 0.00238, 0.0, 0.00071, 0.0, 0.0007, 0.0, 0.00064, 0.0],
        99.999: [ 2.015, 1.29439, 0.45693, 0.70056, 0.45515, 0.62367, 0.45665, 0.574, 0.44846, 0.55042, 0.6195],
        "std_dev": [ 0.1644795868976945, 0.11432684690392436, 0.0686085533598881, 0.1012923823259789, 0.06942697787632549, 0.09637347003804957, 0.07178186491873813, 0.08829919238461785, 0.06838269252500273, 0.08489221834344339, 0.08828317104665424],
        "mean": [ 1.200586290931244, 0.7745734926995385, 0.13121265946721775, 0.28610903922048436, 0.13315419576726892, 0.2094184665635165, 0.1379354620042445, 0.18183810096947914, 0.1319454381771936, 0.1711662322649664, 0.11659769222520298],
    },
    30: {
        5: [ 0.96233, 0.60741, 0.02684, 0.13735, 0.02703, 0.0618, 0.0273, 0.04408, 0.02884, 0.03842, 0.02779, 0.03539, 0.02723, 0.03403, 0.02677, 0.00933],
        95: [ 1.426, 0.92711, 0.20942, 0.41527, 0.20886, 0.31861, 0.21161, 0.28917, 0.22348, 0.27425, 0.21433, 0.25526, 0.20979, 0.24952, 0.2072, 0.278],
        0.1: [ 0.75967, 0.4736, 0.00326, 0.03649, 0.00328, 0.00913, 0.00343, 0.00623, 0.00351, 0.00539, 0.00346, 0.00497, 0.00342, 0.00477, 0.00327, 0.00033],
        99.9: [ 1.68833, 1.10054, 0.31509, 0.53834, 0.31116, 0.43773, 0.31498, 0.40951, 0.32536, 0.39424, 0.31734, 0.37043, 0.31061, 0.36383, 0.30768, 0.45067],
        0.01: [ 0.68567, 0.42631, 0.0, 0.01239, 0.0, 0.00291, 0.0, 0.00197, 0.0, 0.0017, 0.0, 0.00157, 0.0, 0.0015, 0.0, 0.0],
        99.99: [ 1.865, 1.20513, 0.36108, 0.5891, 0.3535, 0.48657, 0.35709, 0.45715, 0.36569, 0.44054, 0.36178, 0.41716, 0.35212, 0.40999, 0.34982, 0.52667],
        0.001: [ 0.61467, 0.38301, 0.0, 0.00393, 0.0, 0.00088, 0.0, 0.00062, 0.0, 0.00054, 0.0, 0.0005, 0.0, 0.00048, 0.0, 0.0],
        99.999: [ 2.015, 1.28514, 0.39977, 0.63166, 0.38817, 0.52661, 0.39148, 0.49495, 0.39775, 0.47716, 0.41, 0.45572, 0.38663, 0.44715, 0.38728, 0.6],
        "std_dev": [ 0.14007808282676615, 0.09658145656885171, 0.05620793298865717, 0.08436481865030923, 0.05592619698969352, 0.07775917771108906, 0.05666550160263995, 0.07472807189718422, 0.05976412494228826, 0.07216318446862328, 0.0573321177540202, 0.06737795047039233, 0.05611136761450671, 0.06608226024276644, 0.05547736162985284, 0.08504655079619773],
        "mean": [ 1.2002046476598136, 0.7701756915293981, 0.10672603302239946, 0.2736331105149904, 0.10687697945712851, 0.18493149427084293, 0.10806175486960326, 0.15684939164168468, 0.11443254462957282, 0.1443968332559008, 0.10973297981816704, 0.1339586746313684, 0.10744542239834769, 0.13010908586627842, 0.10593214051444012, 0.11536600614403454],
    },
    40: {
        5: [ 0.985, 0.62111, 0.02311, 0.14885, 0.0235, 0.06486, 0.02337, 0.04227, 0.0236, 0.03498, 0.02493, 0.03231, 0.02477, 0.02989, 0.02375, 0.02861, 0.02354, 0.02804, 0.02329, 0.02736, 0.00625],
        95: [ 1.40075, 0.90677, 0.18135, 0.39088, 0.18179, 0.29322, 0.18058, 0.25347, 0.18294, 0.23882, 0.19231, 0.23472, 0.1918, 0.22034, 0.18369, 0.21159, 0.18169, 0.20922, 0.18049, 0.20535, 0.20225],
        0.1: [ 0.796, 0.49532, 0.00274, 0.05267, 0.00287, 0.01009, 0.00277, 0.00602, 0.00297, 0.00493, 0.003, 0.00452, 0.00303, 0.00418, 0.00282, 0.00399, 0.00295, 0.00392, 0.00276, 0.00382, 0.0],
        99.9: [ 1.65125, 1.07325, 0.27877, 0.50024, 0.27304, 0.39896, 0.27174, 0.35712, 0.27547, 0.34385, 0.285, 0.34028, 0.28441, 0.32354, 0.27523, 0.31191, 0.27236, 0.30993, 0.27064, 0.30471, 0.34175],
        0.01: [ 0.72325, 0.44907, 0.0, 0.02091, 0.0, 0.00321, 0.0, 0.00191, 0.0, 0.00156, 0.0, 0.00143, 0.0, 0.00133, 0.0, 0.00126, 0.0, 0.00125, 0.0, 0.00121, 0.0],
        99.99: [ 1.841, 1.18663, 0.32375, 0.54689, 0.31229, 0.44367, 0.31141, 0.40104, 0.31575, 0.38737, 0.3242, 0.38256, 0.32242, 0.36776, 0.31496, 0.35426, 0.3117, 0.35297, 0.30905, 0.34766, 0.409],
        0.001: [ 0.615, 0.39192, 0.0, 0.00689, 0.0, 0.00101, 0.0, 0.00061, 0.0, 0.00051, 0.0, 0.00045, 0.0, 0.00042, 0.0, 0.0004, 0.0, 0.00039, 0.0, 0.00039, 0.0],
        99.999: [ 2.00975, 1.28411, 0.36394, 0.5874, 0.34627, 0.481, 0.34599, 0.43869, 0.35013, 0.42325, 0.35902, 0.41775, 0.35335, 0.41088, 0.35187, 0.39052, 0.34779, 0.39004, 0.34229, 0.38619, 0.5015],
        "std_dev": [ 0.1256838042596897, 0.08624029973388889, 0.0488031150869123, 0.073521519832258, 0.04871048057538988, 0.06905884268686814, 0.048382151835161126, 0.06422365153988988, 0.049046362323977445, 0.06230558372140294, 0.05145357265717392, 0.06199487457918859, 0.05134380133111481, 0.05844292683858891, 0.049211434284259746, 0.0561778788432159, 0.04865789270023227, 0.05565099636289747, 0.048363513728984946, 0.05469011007366725, 0.062282380603002736],
        "mean": [ 1.2000044500722438, 0.7682702815698738, 0.09203864086474557, 0.2680917418077898, 0.0928899641548034, 0.17561241688871787, 0.09228391456236364, 0.141405429251508, 0.09329992631600173, 0.1278221966631171, 0.09845023078996198, 0.12268466004303846, 0.09801533310162427, 0.11466499579594534, 0.0938795247492575, 0.10992574012810165, 0.09287553947130722, 0.10827569746285808, 0.09216742930418022, 0.10605803590188662, 0.08188400769635228],
    },
    50: {
        5: [ 1.0006, 0.63076, 0.02058, 0.15777, 0.02102, 0.06872, 0.02094, 0.04221, 0.02089, 0.03309, 0.02107, 0.0294, 0.02206, 0.02803, 0.02252, 0.02664, 0.02165, 0.02532, 0.02116, 0.02457, 0.02104, 0.02429, 0.021, 0.02382, 0.02053, 0.0066],
        95: [ 1.3842, 0.89365, 0.16329, 0.37428, 0.16251, 0.27711, 0.16205, 0.23588, 0.16137, 0.21476, 0.16336, 0.20636, 0.17051, 0.20539, 0.17526, 0.19843, 0.16774, 0.18882, 0.16347, 0.18405, 0.16241, 0.18309, 0.16234, 0.18014, 0.15974, 0.2026],
        0.1: [ 0.8246, 0.51204, 0.00268, 0.06627, 0.00275, 0.01157, 0.00274, 0.00609, 0.00274, 0.00467, 0.00279, 0.00413, 0.00289, 0.00392, 0.00295, 0.00372, 0.00283, 0.00353, 0.00277, 0.00343, 0.00279, 0.0034, 0.00275, 0.00333, 0.00268, 0.0002],
        99.9: [ 1.6262, 1.05501, 0.25624, 0.47331, 0.24563, 0.37299, 0.24479, 0.33037, 0.2441, 0.3077, 0.24759, 0.30096, 0.25585, 0.30037, 0.26049, 0.29289, 0.2528, 0.28154, 0.24617, 0.2738, 0.24574, 0.27364, 0.24512, 0.26937, 0.24209, 0.3376],
        0.01: [ 0.7516, 0.46605, 0.0, 0.03074, 0.0, 0.00368, 0.0, 0.00193, 0.0, 0.00148, 0.00022, 0.0013, 0.0, 0.00124, 0.0, 0.00117, 0.0, 0.00112, 0.0, 0.00109, 0.0002, 0.00107, 0.0, 0.00105, 0.0, 0.0],
        99.99: [ 1.8244, 1.1746, 0.30125, 0.51681, 0.28235, 0.41397, 0.28138, 0.37149, 0.2813, 0.34807, 0.28614, 0.34151, 0.2928, 0.33984, 0.29619, 0.33195, 0.29045, 0.324, 0.28209, 0.31266, 0.28417, 0.31308, 0.28124, 0.3078, 0.28108, 0.4032],
        0.001: [ 0.622, 0.39231, 0.0, 0.011, 0.0, 0.00117, 0.0, 0.00062, 0.0, 0.00047, 0.0, 0.00042, 0.0, 0.0004, 0.0, 0.00037, 0.0, 0.00035, 0.0, 0.00035, 0.0, 0.00034, 0.0, 0.00034, 0.0, 0.0],
        99.999: [ 2.0006, 1.27916, 0.34126, 0.5556, 0.31437, 0.44978, 0.31331, 0.40858, 0.31593, 0.38347, 0.32667, 0.37593, 0.32524, 0.3743, 0.3275, 0.36481, 0.3283, 0.37131, 0.31307, 0.3455, 0.32313, 0.34726, 0.31282, 0.34061, 0.33144, 0.5],
        "std_dev": [ 0.11585601704428589, 0.07931039169042939, 0.04410850381890864, 0.06580328542278098, 0.04356113272205226, 0.06302799654646377, 0.04344277543862622, 0.0587654190113073, 0.04325207086016376, 0.055432041416630065, 0.04381873487812957, 0.054193575147669294, 0.04567795533493363, 0.05439090032473639, 0.046958913292407214, 0.052736933393605, 0.04497150258862254, 0.050250394543536366, 0.04380537102943327, 0.049010220224331144, 0.04353878732923226, 0.04882517346710692, 0.04351178854953644, 0.048066063746765635, 0.042882282728614135, 0.062235457435888056],
        "mean": [ 1.199906420635955, 0.7672156712192271, 0.08237029499604381, 0.26496944430437913, 0.08300349767276544, 0.17053512211772306, 0.08270048722673563, 0.1341858457357938, 0.0823911970138467, 0.11697415276854803, 0.08326914658111444, 0.109353674473439, 0.08712663475635604, 0.1070787605392118, 0.08930135906165873, 0.10272191877189832, 0.08558878465913854, 0.09770154670651601, 0.08347541514040187, 0.09508631163836086, 0.08296368825656311, 0.09432376279283368, 0.08284923334701592, 0.09266562979472322, 0.08134155562241362, 0.08294435576958928],
    },
    60: {
        5: [ 1.0125, 0.63852, 0.01874, 0.16428, 0.019, 0.07258, 0.01915, 0.04301, 0.01908, 0.03235, 0.01904, 0.02775, 0.01922, 0.02562, 0.01995, 0.02488, 0.02062, 0.0242, 0.02023, 0.02301, 0.01955, 0.02226, 0.01922, 0.02182, 0.01917, 0.02168, 0.01917, 0.02139, 0.01879, 0.02104, 0.00517],
        95: [ 1.37183, 0.88416, 0.15074, 0.36246, 0.14746, 0.26519, 0.14853, 0.22382, 0.14754, 0.20151, 0.1473, 0.18865, 0.14897, 0.18342, 0.15454, 0.18351, 0.1599, 0.18209, 0.15699, 0.17321, 0.1515, 0.16778, 0.14879, 0.16478, 0.14824, 0.16449, 0.14877, 0.16251, 0.14613, 0.16061, 0.16483],
        0.1: [ 0.84667, 0.52688, 0.00231, 0.07551, 0.0024, 0.01354, 0.00241, 0.00631, 0.00241, 0.0046, 0.00236, 0.00389, 0.00249, 0.00359, 0.00248, 0.0035, 0.0026, 0.00339, 0.00255, 0.00321, 0.0025, 0.00311, 0.00239, 0.00305, 0.00249, 0.00303, 0.00238, 0.00299, 0.00236, 0.00294, 0.00017],
        99.9: [ 1.60633, 1.04079, 0.24121, 0.45377, 0.22399, 0.35372, 0.22504, 0.31063, 0.22416, 0.28717, 0.22348, 0.27344, 0.22675, 0.26985, 0.23329, 0.27106, 0.239, 0.27057, 0.23597, 0.25879, 0.23047, 0.25116, 0.22511, 0.24656, 0.22585, 0.24701, 0.22611, 0.24417, 0.22177, 0.24237, 0.28133],
        0.01: [ 0.771, 0.48104, 0.0, 0.03836, 0.0, 0.00432, 0.0, 0.002, 0.0, 0.00145, 0.0, 0.00124, 0.0, 0.00113, 0.0, 0.00107, 0.0, 0.00107, 0.0, 0.00102, 0.0, 0.00098, 0.0, 0.00096, 0.0, 0.00095, 0.0, 0.00094, 0.0, 0.00093, 0.0],
        99.99: [ 1.81017, 1.16469, 0.28549, 0.49501, 0.25823, 0.39227, 0.25899, 0.34806, 0.25911, 0.32514, 0.25879, 0.31092, 0.26392, 0.30778, 0.26786, 0.31001, 0.27262, 0.30812, 0.26969, 0.29599, 0.27023, 0.28764, 0.25898, 0.28244, 0.26467, 0.28329, 0.26201, 0.28016, 0.25521, 0.28075, 0.34333],
        0.001: [ 0.63067, 0.39753, 0.0, 0.01484, 0.0, 0.00137, 0.0, 0.00064, 0.0, 0.00045, 0.0, 0.00039, 0.0, 0.00036, 0.0, 0.00037, 0.0, 0.00034, 0.0, 0.00033, 0.0, 0.00031, 0.0, 0.0003, 0.0, 0.0003, 0.0, 0.0003, 0.0, 0.0003, 0.0],
        99.999: [ 1.99667, 1.27509, 0.32503, 0.53223, 0.2883, 0.42566, 0.28872, 0.38098, 0.29181, 0.3599, 0.2953, 0.34574, 0.31107, 0.34124, 0.29749, 0.34817, 0.30165, 0.34093, 0.29836, 0.32726, 0.33226, 0.31929, 0.2901, 0.31458, 0.30823, 0.31532, 0.29957, 0.31184, 0.28478, 0.32917, 0.4745],
        "std_dev": [ 0.10847732931912153, 0.07404236898542912, 0.04089470357246045, 0.060279663495733445, 0.03956622203144364, 0.058315426277901615, 0.03984384170492117, 0.05476289763764866, 0.039568750741002946, 0.05153314440771936, 0.039499641902402884, 0.04921185104678147, 0.03998081804800849, 0.04839046467344727, 0.041436086870757, 0.048701527835043286, 0.042839865771181984, 0.04850957349307759, 0.042083456624326394, 0.04617153678819423, 0.04065221785240769, 0.04474256875642022, 0.039898683347052384, 0.04396009343784075, 0.03977995631517226, 0.043931624227812736, 0.039926111886939494, 0.043416791103570045, 0.03922033513695396, 0.04295946751626444, 0.050851399516256446],
        "mean": [ 1.1997950331720049, 0.7665638012651648, 0.07548542282972927, 0.26288479047056396, 0.07517965229320521, 0.16716812645798992, 0.07571635407511318, 0.12963364475634317, 0.07529499161704854, 0.11118198956698976, 0.07514921829180606, 0.1011698106717842, 0.07591242350663814, 0.09647433867227143, 0.07884323791017615, 0.09540925919944998, 0.08160836333092923, 0.0937996362067251, 0.0800666734957316, 0.08927608901622582, 0.0772611797234559, 0.08644024045530574, 0.07593445831436259, 0.08485177609148911, 0.0756737370618472, 0.08450687620788569, 0.07582743710371255, 0.0834357969672141, 0.07443383428704656, 0.0823286552741688, 0.06676791064282225],
    },
    70: {
        5: [ 1.02114, 0.64475, 0.01728, 0.16884, 0.01742, 0.07655, 0.01777, 0.04418, 0.01767, 0.03215, 0.01765, 0.0269, 0.01762, 0.02423, 0.01777, 0.02291, 0.01833, 0.02247, 0.01899, 0.02215, 0.01899, 0.02138, 0.01837, 0.02056, 0.01798, 0.02008, 0.01778, 0.01982, 0.01774, 0.01974, 0.01778, 0.01955, 0.01752, 0.01923, 0.0173, 0.00529],
        95: [ 1.362, 0.87674, 0.14103, 0.35387, 0.13535, 0.2555, 0.13765, 0.2146, 0.13711, 0.19186, 0.1364, 0.17814, 0.13628, 0.16964, 0.13774, 0.1663, 0.14221, 0.16682, 0.14723, 0.16737, 0.14869, 0.162, 0.14257, 0.156, 0.13926, 0.15234, 0.13746, 0.15041, 0.13721, 0.15035, 0.13804, 0.14923, 0.13593, 0.147, 0.13481, 0.16614],
        0.1: [ 0.86229, 0.53843, 0.00227, 0.08172, 0.0023, 0.0161, 0.00235, 0.00664, 0.00234, 0.00459, 0.00237, 0.00378, 0.00233, 0.00339, 0.00237, 0.00321, 0.00242, 0.00314, 0.00251, 0.0031, 0.00254, 0.00299, 0.00243, 0.00288, 0.00238, 0.00281, 0.00236, 0.00277, 0.00237, 0.00277, 0.00238, 0.00273, 0.00232, 0.00269, 0.00229, 0.00014],
        99.9: [ 1.59086, 1.02951, 0.22957, 0.43918, 0.20673, 0.33782, 0.20928, 0.29566, 0.20842, 0.27176, 0.208, 0.25712, 0.2073, 0.24814, 0.21043, 0.24646, 0.21588, 0.24806, 0.22191, 0.24952, 0.22464, 0.24269, 0.2165, 0.2356, 0.21206, 0.22909, 0.2087, 0.22601, 0.21034, 0.22641, 0.21085, 0.22515, 0.20694, 0.22226, 0.20636, 0.28029],
        0.01: [ 0.78686, 0.49228, 0.00021, 0.04439, 0.00021, 0.00519, 0.00022, 0.00211, 0.00022, 0.00146, 0.00039, 0.0012, 0.00022, 0.00107, 0.00033, 0.00101, 0.00023, 0.00099, 0.00024, 0.00097, 0.00042, 0.00095, 0.00023, 0.0009, 0.00022, 0.00089, 0.00022, 0.00087, 0.00033, 0.00088, 0.0004, 0.00087, 0.00022, 0.00085, 0.00021, 0.0],
        99.99: [ 1.79643, 1.15514, 0.27261, 0.47885, 0.23949, 0.37364, 0.2416, 0.33092, 0.24041, 0.30684, 0.24265, 0.29211, 0.24067, 0.2831, 0.24615, 0.2823, 0.24854, 0.28409, 0.25526, 0.28456, 0.25818, 0.27745, 0.24869, 0.27428, 0.24664, 0.26318, 0.24142, 0.2598, 0.25075, 0.26008, 0.24615, 0.25893, 0.23911, 0.25573, 0.24317, 0.341],
        0.001: [ 0.641, 0.40177, 0.0, 0.01849, 0.0, 0.00165, 0.0, 0.00067, 0.0, 0.00046, 0.0, 0.00038, 0.0, 0.00034, 0.0, 0.00032, 0.0, 0.00031, 0.0, 0.00031, 0.0, 0.0003, 0.0, 0.00028, 0.0, 0.00028, 0.0, 0.00028, 0.0, 0.00028, 0.0, 0.00027, 0.0, 0.00027, 0.0, 0.0],
        99.999: [ 1.98929, 1.27075, 0.30988, 0.51603, 0.2683, 0.40521, 0.27029, 0.36221, 0.26995, 0.3385, 0.28285, 0.32426, 0.27786, 0.31698, 0.29961, 0.31417, 0.27882, 0.31819, 0.28847, 0.31489, 0.29158, 0.30759, 0.27519, 0.32229, 0.28524, 0.29384, 0.27454, 0.2898, 0.30028, 0.28993, 0.29071, 0.28878, 0.26864, 0.2858, 0.30753, 0.45943],
        "std_dev": [ 0.10285257001987054, 0.06990518052291728, 0.038412647666526266, 0.05631237239012467, 0.0363400250439689, 0.054228977098579025, 0.036929518912335944, 0.05156275447152133, 0.03679226251532327, 0.04857093101726372, 0.03659508761839197, 0.04620334667145834, 0.03655975738579279, 0.04454900465004169, 0.03698667521798811, 0.04401545001819817, 0.03816115194978453, 0.044343782281541004, 0.03947698734035504, 0.04462959257783958, 0.03993611662357533, 0.04324156797852145, 0.03825103637577571, 0.04169462791791582, 0.03737009371266489, 0.04068961618537642, 0.0368614039198378, 0.04017313092030677, 0.03684619670263958, 0.040191503873507735, 0.03706810485386983, 0.03990629236486346, 0.036483122134927964, 0.039334522167483205, 0.0362342243388592, 0.05114558622712205],
        "mean": [ 1.1995230812332602, 0.7660306053552375, 0.07008864220093694, 0.2613476165253962, 0.0689233275887079, 0.16468888228191322, 0.07021564274376994, 0.1264058161179121, 0.06986763514320818, 0.10716777514453152, 0.06959125421399097, 0.09637526836806078, 0.06952111158522004, 0.09001315241115848, 0.07019523056956098, 0.08698163383019861, 0.07245978425540466, 0.08648774864111512, 0.07510781721602897, 0.08609083922051011, 0.07545496392915622, 0.08325846736699542, 0.0726120545006893, 0.0801458181038083, 0.07102036486991753, 0.0782644245438234, 0.07014554012563631, 0.0772697074320454, 0.07001309165205902, 0.07710760237954825, 0.0702925448726329, 0.07646167384870911, 0.06924279042602892, 0.07527842292370372, 0.0685506672668014, 0.06772230952979019],
    },
    80: {
        5: [ 1.02737, 0.64944, 0.01609, 0.17235, 0.01619, 0.08037, 0.01657, 0.04552, 0.01656, 0.03229, 0.0165, 0.02641, 0.01647, 0.0234, 0.01645, 0.0217, 0.01661, 0.02083, 0.01703, 0.02058, 0.01764, 0.02045, 0.01781, 0.02002, 0.0175, 0.01929, 0.01698, 0.01879, 0.01672, 0.01844, 0.01661, 0.01827, 0.01658, 0.01822, 0.01662, 0.01812, 0.01643, 0.01784, 0.01619, 0.01769, 0.0045],
        95: [ 1.35438, 0.87097, 0.13288, 0.34738, 0.12608, 0.2473, 0.12854, 0.20713, 0.12867, 0.18442, 0.1279, 0.17026, 0.12757, 0.16102, 0.12746, 0.15508, 0.12876, 0.15292, 0.13241, 0.15354, 0.13695, 0.15457, 0.13977, 0.15327, 0.13611, 0.14683, 0.13196, 0.143, 0.12966, 0.14046, 0.12845, 0.13916, 0.12836, 0.13922, 0.12922, 0.13867, 0.12762, 0.13672, 0.12605, 0.1359, 0.14263],
        0.1: [ 0.87475, 0.54784, 0.00206, 0.0864, 0.00209, 0.01939, 0.00213, 0.00704, 0.00217, 0.00464, 0.00214, 0.00372, 0.00213, 0.00329, 0.00212, 0.00304, 0.00219, 0.00292, 0.00218, 0.00288, 0.0023, 0.00286, 0.00229, 0.0028, 0.00229, 0.00269, 0.00218, 0.00262, 0.00216, 0.00258, 0.00215, 0.00255, 0.00219, 0.00255, 0.00214, 0.00253, 0.00213, 0.0025, 0.00208, 0.00247, 0.00012],
        99.9: [ 1.57837, 1.02041, 0.21969, 0.42902, 0.19328, 0.32447, 0.19546, 0.28362, 0.19574, 0.25977, 0.19486, 0.24508, 0.19467, 0.23469, 0.19409, 0.22842, 0.19718, 0.22789, 0.20191, 0.2291, 0.20809, 0.23117, 0.21172, 0.23145, 0.20679, 0.22117, 0.20396, 0.21674, 0.19731, 0.21191, 0.19541, 0.20997, 0.19774, 0.21041, 0.19783, 0.21031, 0.19454, 0.20734, 0.1925, 0.2072, 0.245],
        0.01: [ 0.79838, 0.50172, 0.0, 0.04966, 0.0, 0.00646, 0.0, 0.00222, 0.0, 0.00147, 0.0, 0.00118, 0.0, 0.00104, 0.0, 0.00096, 0.0, 0.00092, 0.0, 0.00091, 0.0, 0.0009, 0.0, 0.00088, 0.0, 0.00085, 0.0, 0.00083, 0.0, 0.00081, 0.0, 0.00081, 0.0, 0.00081, 0.0, 0.0008, 0.0, 0.00079, 0.0, 0.00078, 0.0],
        99.99: [ 1.782, 1.1452, 0.26138, 0.46799, 0.22463, 0.35849, 0.22591, 0.31736, 0.22601, 0.29288, 0.22564, 0.27933, 0.22618, 0.26814, 0.22573, 0.26171, 0.23178, 0.26202, 0.23317, 0.26234, 0.24239, 0.26416, 0.24337, 0.26569, 0.23797, 0.25381, 0.242, 0.25382, 0.22843, 0.24407, 0.22686, 0.24204, 0.24054, 0.24234, 0.23182, 0.24304, 0.22538, 0.23913, 0.22258, 0.24318, 0.3045],
        0.001: [ 0.64688, 0.40754, 0.0, 0.02291, 0.0, 0.00207, 0.0, 0.0007, 0.0, 0.00047, 0.0, 0.00038, 0.0, 0.00033, 0.0, 0.0003, 0.0, 0.00029, 0.0, 0.00029, 0.0, 0.00029, 0.0, 0.00028, 0.0, 0.00027, 0.0, 0.00026, 0.0, 0.00026, 0.0, 0.00026, 0.0, 0.00025, 0.0, 0.00025, 0.0, 0.00025, 0.0, 0.00025, 0.0],
        99.999: [ 1.984, 1.26791, 0.29738, 0.50412, 0.25341, 0.38788, 0.25347, 0.34716, 0.25469, 0.32297, 0.25647, 0.31339, 0.25943, 0.3024, 0.26069, 0.2939, 0.28976, 0.29435, 0.26398, 0.29121, 0.29001, 0.29279, 0.27206, 0.29919, 0.26574, 0.28095, 0.27955, 0.30722, 0.25651, 0.27397, 0.26167, 0.27119, 0.29362, 0.27097, 0.27796, 0.27462, 0.25469, 0.26802, 0.25027, 0.30617, 0.44213],
        "std_dev": [ 0.09853319697380297, 0.06669150445980303, 0.03632692832789364, 0.05329945460918929, 0.03388148278191862, 0.05063039187962177, 0.03449784383740368, 0.048881511118356015, 0.034538753434888506, 0.04619511061597051, 0.03432413023641067, 0.043902816593296756, 0.0342391140504576, 0.042123809607617024, 0.034205866790451435, 0.04090671486994288, 0.03457906281772196, 0.04057479324486056, 0.03555578168947146, 0.04086237527158287, 0.03676498280454805, 0.041232395740398026, 0.037563212042439315, 0.041012133079292586, 0.03654205497170288, 0.03924111310858344, 0.0354887321136906, 0.03825412693615443, 0.03479366781751557, 0.03755160622768339, 0.034455846818660256, 0.03720837742848884, 0.0344923522483952, 0.03724601118821303, 0.03471707375773967, 0.03712359277465369, 0.034260978533012486, 0.03660526338408506, 0.033854020941270435, 0.03642704350345188, 0.04404405010655183],
        "mean": [ 1.1992301845464959, 0.7655794627764031, 0.06560252222043077, 0.2602587028125359, 0.06412002206169064, 0.16278909239790876, 0.06552702659096238, 0.12398399416342973, 0.06552057719754592, 0.10422596003723478, 0.0651760409740162, 0.09292916133778494, 0.06502874574015079, 0.08597955559271177, 0.06499731573172837, 0.08165738602999197, 0.06559881230447934, 0.07963703040234149, 0.06739853789187801, 0.07942071002802605, 0.0697913768983556, 0.07951185745361361, 0.07091491551865217, 0.07837466230327642, 0.06925805025897104, 0.07528895249785877, 0.06719491300824103, 0.07336200719615797, 0.06607151173032047, 0.07204532431892488, 0.06552280242916066, 0.07137234984679144, 0.06545462342964484, 0.07130036962269339, 0.06575886611113879, 0.07095773803744977, 0.06499354306103855, 0.06991630863738453, 0.06412792166390803, 0.06944055644842252, 0.05778131768353812],
    },
    90: {
        5: [ 1.03189, 0.65308, 0.01509, 0.17548, 0.01525, 0.08371, 0.01551, 0.04708, 0.01567, 0.03259, 0.01559, 0.02616, 0.01555, 0.02285, 0.01551, 0.02096, 0.01551, 0.01981, 0.01565, 0.01921, 0.016, 0.01904, 0.01654, 0.01902, 0.01679, 0.01882, 0.01666, 0.0183, 0.01621, 0.01774, 0.01592, 0.01739, 0.01569, 0.01717, 0.01564, 0.01703, 0.01563, 0.017, 0.01567, 0.01694, 0.01556, 0.01672, 0.01535, 0.01651, 0.01524, 0.00456],
        95: [ 1.348, 0.86621, 0.12572, 0.34209, 0.11924, 0.24054, 0.12052, 0.20078, 0.12149, 0.17829, 0.12099, 0.16394, 0.12034, 0.1545, 0.12019, 0.14786, 0.1202, 0.14354, 0.12132, 0.14214, 0.12442, 0.14272, 0.12836, 0.14401, 0.13115, 0.14451, 0.13114, 0.13993, 0.1262, 0.13556, 0.12363, 0.13282, 0.12184, 0.13097, 0.12103, 0.13004, 0.12102, 0.1302, 0.12184, 0.13007, 0.12067, 0.1283, 0.11918, 0.12703, 0.1188, 0.14389],
        0.1: [ 0.88456, 0.55547, 0.002, 0.09079, 0.00202, 0.02308, 0.00207, 0.00758, 0.00208, 0.00473, 0.00209, 0.00371, 0.00208, 0.00322, 0.00207, 0.00292, 0.00207, 0.00277, 0.00211, 0.00269, 0.00215, 0.00266, 0.0022, 0.00266, 0.00225, 0.00263, 0.00222, 0.00255, 0.00216, 0.00247, 0.00214, 0.00243, 0.00209, 0.00239, 0.00208, 0.00237, 0.0021, 0.00237, 0.00209, 0.00237, 0.00209, 0.00233, 0.00206, 0.0023, 0.00203, 0.00011],
        99.9: [ 1.56811, 1.01294, 0.21018, 0.42151, 0.18352, 0.31345, 0.18374, 0.27315, 0.18515, 0.24989, 0.18446, 0.23462, 0.18366, 0.2248, 0.18356, 0.21718, 0.18333, 0.21268, 0.18605, 0.21291, 0.19031, 0.2137, 0.19564, 0.21654, 0.19907, 0.21884, 0.20046, 0.21185, 0.19185, 0.2081, 0.18968, 0.20143, 0.18585, 0.19817, 0.18449, 0.19698, 0.18725, 0.19743, 0.18645, 0.19854, 0.18415, 0.19508, 0.18249, 0.19323, 0.18314, 0.24511],
        0.01: [ 0.80722, 0.50804, 0.00028, 0.05485, 0.00029, 0.00792, 0.00035, 0.0024, 0.0003, 0.00149, 0.00036, 0.00117, 0.00035, 0.00102, 0.00029, 0.00091, 0.00029, 0.00088, 0.0004, 0.00084, 0.00037, 0.00084, 0.00032, 0.00084, 0.00039, 0.00083, 0.00031, 0.00081, 0.00031, 0.00078, 0.0004, 0.00077, 0.00029, 0.00076, 0.00029, 0.00075, 0.00041, 0.00075, 0.0003, 0.00075, 0.00036, 0.00074, 0.00035, 0.00073, 0.00028, 0.0],
        99.99: [ 1.769, 1.13684, 0.25018, 0.46012, 0.21347, 0.346, 0.2125, 0.30495, 0.21437, 0.28152, 0.2134, 0.26609, 0.21326, 0.25772, 0.21304, 0.24951, 0.21314, 0.24408, 0.21984, 0.2455, 0.22064, 0.24517, 0.22757, 0.2497, 0.22909, 0.25157, 0.23161, 0.24379, 0.22112, 0.24799, 0.22634, 0.23304, 0.21543, 0.22907, 0.21469, 0.22764, 0.2329, 0.22807, 0.21767, 0.23158, 0.21377, 0.22576, 0.2117, 0.22343, 0.21882, 0.30322],
        0.001: [ 0.65367, 0.4121, 0.0, 0.0279, 0.0, 0.00257, 0.0, 0.00074, 0.0, 0.00047, 0.0, 0.00038, 0.0, 0.00032, 0.0, 0.00029, 0.0, 0.00028, 0.0, 0.00027, 0.0, 0.00027, 0.0, 0.00027, 0.0, 0.00026, 0.0, 0.00025, 0.0, 0.00024, 0.0, 0.00025, 0.0, 0.00024, 0.0, 0.00024, 0.0, 0.00024, 0.0, 0.00024, 0.0, 0.00023, 0.0, 0.00023, 0.0, 0.0],
        99.999: [ 1.979, 1.26437, 0.2861, 0.49694, 0.24153, 0.37583, 0.23846, 0.33337, 0.24129, 0.31046, 0.24181, 0.29541, 0.24223, 0.2961, 0.24107, 0.28763, 0.24615, 0.27673, 0.28414, 0.27599, 0.25117, 0.27388, 0.26687, 0.28394, 0.25589, 0.28092, 0.2656, 0.27303, 0.2474, 0.28409, 0.2968, 0.26327, 0.24322, 0.25898, 0.24838, 0.25539, 0.28752, 0.25662, 0.25698, 0.27059, 0.24389, 0.25513, 0.23825, 0.25123, 0.29837, 0.43944],
        "std_dev": [ 0.09498372909753593, 0.06408909436228305, 0.034458620470906476, 0.0507764860944925, 0.032076402949056795, 0.04759844429311936, 0.0323597393589157, 0.046481231589917044, 0.03260411606423302, 0.04418454894544713, 0.03248092749297168, 0.04198476365628616, 0.0322920683793858, 0.0402732039957564, 0.03225818562608094, 0.038896595239352016, 0.03226342914184631, 0.037985196915230106, 0.03258522720251871, 0.03778727758127818, 0.03341752501443082, 0.03802581483219264, 0.03446375442226942, 0.038450238550029756, 0.035226099721668305, 0.038694457845271577, 0.03529287073284626, 0.03744601269099978, 0.03388177955831244, 0.03632803163601302, 0.03322399862664645, 0.035544689802498836, 0.032714143857988806, 0.03503394959500662, 0.03247486156619183, 0.034798318791130746, 0.03253560794345095, 0.034860571022503034, 0.03272946514470902, 0.03486030092294054, 0.03238657807156906, 0.03437156535497627, 0.03200062248436463, 0.03404544405165921, 0.03195891605911791, 0.04437330849713467],
        "mean": [ 1.1988653404109817, 0.7651470532013851, 0.061754557650871754, 0.25943477205045784, 0.060539171583137145, 0.16130955702184585, 0.06135108207703439, 0.12205837981935952, 0.061915915448674745, 0.10191570898750835, 0.06162076701855885, 0.09023463185086968, 0.061338836331387925, 0.08300737260454913, 0.061228755760396505, 0.07823553194819144, 0.0612723190816112, 0.07516306273579328, 0.061799003497139834, 0.07377559338191683, 0.0632912028649296, 0.07366578660537448, 0.06538731540652329, 0.07402967241496827, 0.06662787955124981, 0.07378016738491437, 0.0663518340113519, 0.07161569630060817, 0.06419088517567527, 0.06940694901948878, 0.06294624176386734, 0.06803421613557788, 0.06206722259518999, 0.06711980345398023, 0.06171891600778187, 0.06662578414874099, 0.06168660892270912, 0.06660198247946582, 0.06199295216673278, 0.06645906336012196, 0.06146441798130132, 0.0655632715878444, 0.060650444574308404, 0.0648767337142286, 0.06037337016380196, 0.05854479475770266],
    },
    100: {
        5: [ 1.0351, 0.65587, 0.01423, 0.17838, 0.01452, 0.08649, 0.01457, 0.04882, 0.01487, 0.03302, 0.0148, 0.02606, 0.01475, 0.02248, 0.01474, 0.02042, 0.01471, 0.01914, 0.01472, 0.01832, 0.01485, 0.01789, 0.01514, 0.01777, 0.01559, 0.01781, 0.01591, 0.01772, 0.01591, 0.01741, 0.01562, 0.01691, 0.01524, 0.01653, 0.01501, 0.01627, 0.01489, 0.01612, 0.01482, 0.01601, 0.01482, 0.01598, 0.01487, 0.01597, 0.01478, 0.01578, 0.01458, 0.0156, 0.01442, 0.01554, 0.004],
        95: [ 1.3426, 0.86222, 0.11938, 0.33751, 0.11407, 0.23522, 0.11333, 0.19521, 0.11517, 0.17309, 0.11514, 0.15879, 0.11452, 0.14891, 0.11415, 0.14232, 0.11395, 0.13737, 0.11408, 0.13416, 0.11502, 0.13323, 0.11767, 0.13371, 0.12109, 0.13524, 0.12365, 0.13602, 0.12545, 0.13445, 0.1218, 0.12945, 0.11853, 0.12677, 0.11669, 0.12454, 0.11529, 0.12322, 0.11475, 0.1225, 0.11478, 0.1227, 0.11552, 0.12281, 0.11475, 0.12123, 0.11336, 0.11996, 0.11244, 0.11974, 0.1275],
        0.1: [ 0.8937, 0.56169, 0.00188, 0.09537, 0.00193, 0.0264, 0.00193, 0.00821, 0.00197, 0.00484, 0.00198, 0.00371, 0.00196, 0.00316, 0.00196, 0.00286, 0.00196, 0.00269, 0.00195, 0.00256, 0.002, 0.0025, 0.00201, 0.00248, 0.00208, 0.0025, 0.00211, 0.00247, 0.00212, 0.00243, 0.00208, 0.00236, 0.00203, 0.00231, 0.00199, 0.00227, 0.00198, 0.00225, 0.00197, 0.00223, 0.002, 0.00223, 0.00197, 0.00223, 0.00196, 0.0022, 0.00193, 0.00218, 0.00191, 0.00217, 0.0001],
        99.9: [ 1.5586, 1.00642, 0.20084, 0.41534, 0.17633, 0.30483, 0.17341, 0.26407, 0.17569, 0.2413, 0.17571, 0.22633, 0.1749, 0.2155, 0.17485, 0.20871, 0.17411, 0.20292, 0.17424, 0.19972, 0.1766, 0.20025, 0.18085, 0.20086, 0.18432, 0.20438, 0.18814, 0.20621, 0.19245, 0.20526, 0.18628, 0.19583, 0.1857, 0.19362, 0.17897, 0.18912, 0.17609, 0.187, 0.17516, 0.18611, 0.17827, 0.18668, 0.17667, 0.18856, 0.17535, 0.18458, 0.17364, 0.1829, 0.17247, 0.18392, 0.2201],
        0.01: [ 0.8145, 0.5131, 0.0002, 0.06038, 0.00026, 0.00957, 0.00021, 0.00262, 0.00026, 0.00153, 0.0003, 0.00117, 0.00026, 0.00099, 0.00021, 0.00091, 0.00025, 0.00086, 0.00021, 0.00082, 0.00038, 0.00079, 0.00022, 0.00079, 0.00027, 0.0008, 0.00022, 0.00078, 0.00028, 0.00077, 0.00032, 0.00075, 0.00027, 0.00073, 0.00021, 0.00072, 0.00026, 0.00071, 0.00021, 0.0007, 0.00039, 0.00071, 0.00022, 0.0007, 0.00026, 0.00069, 0.00021, 0.00069, 0.00025, 0.00069, 0.0],
        99.99: [ 1.757, 1.12932, 0.23996, 0.45488, 0.20549, 0.33641, 0.20108, 0.29447, 0.20344, 0.27147, 0.20382, 0.25623, 0.20315, 0.24535, 0.20552, 0.23908, 0.2028, 0.2334, 0.20261, 0.2299, 0.20983, 0.23153, 0.21031, 0.23081, 0.21314, 0.2378, 0.21696, 0.23767, 0.22272, 0.23697, 0.21574, 0.22562, 0.22764, 0.22786, 0.21098, 0.21873, 0.20484, 0.21694, 0.20389, 0.21573, 0.2278, 0.21609, 0.20552, 0.22238, 0.20365, 0.21417, 0.20171, 0.21194, 0.20044, 0.21911, 0.2786],
        0.001: [ 0.659, 0.41584, 0.0, 0.03315, 0.0, 0.00305, 0.0, 0.00083, 0.0, 0.00048, 0.0, 0.00037, 0.0, 0.00031, 0.0, 0.00028, 0.0, 0.00026, 0.0, 0.00026, 0.0, 0.00025, 0.0, 0.00025, 0.0, 0.00022, 0.0, 0.00025, 0.0, 0.00025, 0.0, 0.00024, 0.0, 0.00023, 0.0, 0.00023, 0.0, 0.00022, 0.0, 0.00022, 0.0, 0.00023, 0.0, 0.00022, 0.0, 0.00022, 0.0, 0.00021, 0.0, 0.00022, 0.0],
        99.999: [ 1.9765, 1.26232, 0.27638, 0.49264, 0.23205, 0.36491, 0.22652, 0.32202, 0.22972, 0.29897, 0.23213, 0.28451, 0.23153, 0.27406, 0.2511, 0.26937, 0.23506, 0.26986, 0.2356, 0.26226, 0.27763, 0.26172, 0.24088, 0.25934, 0.24065, 0.28476, 0.24276, 0.26591, 0.2533, 0.26862, 0.24428, 0.25213, 0.2597, 0.28417, 0.25264, 0.24649, 0.23354, 0.24745, 0.23552, 0.24382, 0.28184, 0.24388, 0.23781, 0.27167, 0.23304, 0.24411, 0.22976, 0.23932, 0.22844, 0.29493, 0.4253],
        "std_dev": [ 0.09205206182268906, 0.061960488142494385, 0.0327787524044472, 0.04852924754852681, 0.030716914714269426, 0.04516589099208486, 0.0304454557952896, 0.044275963756641126, 0.030905286562552144, 0.04243394259018718, 0.030917434728119805, 0.04040311085558669, 0.030746993277073097, 0.03863373107501763, 0.030649034000591723, 0.03734881076975584, 0.030582004247851086, 0.03627705687869411, 0.030626411207138705, 0.03558293714590452, 0.030897897015121273, 0.0354718156859136, 0.031620559179386415, 0.035663488878274203, 0.03250532571634258, 0.03614766105586389, 0.03319494348991628, 0.036424385937555745, 0.03378015090640423, 0.036065488194352596, 0.03273002612686294, 0.03464347482915386, 0.031917668797452385, 0.033980522775730666, 0.03135579484775615, 0.033344134328251486, 0.030942769415932055, 0.032984326353803425, 0.030799322623475255, 0.03279989456321568, 0.030879141955547914, 0.03287329445650008, 0.031024155802116513, 0.032950330064648696, 0.03080943632178213, 0.03249137939417016, 0.03044858822248168, 0.03215386803421595, 0.030228505554146013, 0.032140410673119006, 0.039403855232722426],
        "mean": [ 1.198464347953969, 0.7647283938534896, 0.05841558997525016, 0.25875700214075537, 0.05777787183223377, 0.1602053085781776, 0.05765924986702958, 0.12046248911450255, 0.0587230234273827, 0.10008521647194923, 0.05857942243111472, 0.0881225585247866, 0.05831240299108356, 0.0805429106024091, 0.058171423282475675, 0.07562999453922517, 0.058025629824132736, 0.07220898707416724, 0.05813407495123435, 0.06994243853255747, 0.058595827460049016, 0.06896679385437275, 0.059859546645054454, 0.06890254286618222, 0.061658283630538754, 0.06944273289205934, 0.06292211769970486, 0.06946197751210315, 0.06340309527285498, 0.06848186561375157, 0.06188230084066179, 0.06622038999735527, 0.060303901589104744, 0.06481243461751103, 0.059399236064223745, 0.06372101845466734, 0.05875654170838824, 0.06308746763679095, 0.05850127642897355, 0.06270499511219811, 0.05849386207431702, 0.0627062999150859, 0.05877739845166965, 0.06268824631434193, 0.05843447488233042, 0.061932519788337026, 0.05768609666019699, 0.0612482018856901, 0.057173382198883724, 0.061092537789193734, 0.05164927098476661],
    }
}


def create_constellation(
        aa_vec: np.ndarray,
        config: DBConfig,
        window=WINDOW_TYPE
        ) -> ConstellationMap:
    """
    The function carries out a windowed fast fourier transformation,
    often called short time FFT (STFT), on the given vector and creates a list
    of frequencies, found by the STFT, and the index of the window they are found

    ...

    Parameter
    ---------
    aa_vec : np.ndarray
        An array of floats representing an amino acid sequence
    window_size : int
        size of a window, related to amino acids -> size 2 means two amino acids
    n_peaks : int, optional
        the number of frequency peaks that are selected from the STFT result,
        if 0 then all are selected (defaults to 0)
    window : str
        The window type for the STFT, read https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.get_window.html#scipy.signal.get_window
        for supported ones (defaults to "boxcar")
    overlap : int
        the overlap between two windows during doing the STFT (defaults to half the window size)

    Returns
    -------
     A ConstellationMap, a list of coordinates, meaning the pairs of
     window index and prominent frequency peak of the window
    """

    if len(aa_vec) < config.window_size:
        aa_vec = np.pad(aa_vec, (0, config.window_size - len(aa_vec)))
        assert len(aa_vec) == config.window_size

    # executing the STFT
    stft_result = signal.stft(
        aa_vec,
        nperseg=config.window_size,
        noverlap=config.overlap,
        window=window
    )

    return stft_to_constellation(*stft_result, config)


def stft_to_constellation(
        frequencies: np.ndarray,
        window_indexes: np.ndarray,
        stft: np.ndarray,
        config: DBConfig
        ) -> ConstellationMap:
    constellation_map: ConstellationMap = []

    # find and collect the most prominent frequencies from STFT per window
    for amplitudes in stft.T:

        # get rid of complex values to make them comparable
        spectrum: np.ndarray = abs(amplitudes)

        peaks: List[Tuple[int, int]] = find_peaks(spectrum, config)

        constellation_map.append(tuple((int(freq_idx), float(spectrum[freq_idx]), quantile) for freq_idx, quantile in peaks))

    return constellation_map


def find_peaks(spectrum: np.ndarray, config: DBConfig) -> List[Tuple[int, int]]:
    lower, upper = sorted((config.significance, 100 - config.significance))
    peaks = []

    for quantile, (tail, op, spectrum_mult) in enumerate(((upper, gt, 1), (lower, lt, -1))):
        quantile_deviations = abs(spectrum - AMPL_QUANTILES[config.window_size][tail]) / AMPL_QUANTILES[config.window_size]["std_dev"]

        if config.selection_method == "none":
            peak_idx = list(range(len(spectrum)))
        elif config.selection_method == "deviation":
            peak_idx, _ = signal.find_peaks(quantile_deviations)
        elif config.selection_method == "absolute":
            peak_idx, _ = signal.find_peaks(spectrum * spectrum_mult)

        tail_idx = np.argwhere(op(spectrum, AMPL_QUANTILES[config.window_size][tail])).flatten()
        selection_idx = np.intersect1d(tail_idx, peak_idx)

        if not len(selection_idx):
            continue

        selected_deviations = quantile_deviations[selection_idx]
        peaks += [*zip(selected_deviations.tolist(), selection_idx.tolist(), [quantile] * len(selected_deviations))]

    peaks = sorted(peaks, reverse=True)
    if config.n_peaks:
        peaks = peaks[:config.n_peaks]

    return [(idx, q) for _, idx, q in peaks if not idx < config.skip_first_k_freqs]
