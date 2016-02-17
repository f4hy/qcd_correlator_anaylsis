def get_heavyq_mass(beta, heavytype):
    if heavytype is None:
        return None
    print beta
    if beta == "4.17":
        #0.44037 0.55046 0.68808 0.86001
        heavymap = {"m0": 0.44037, "m1": 0.55046, "m2": 0.68808, "m3": 0.86001, "m4": float("NAN"), "m5": float("NAN")}
    if beta == "4.35":
        #0.27287 0.34109 0.42636 0.53295 0.66619 0.83273
        heavymap = {"m0": 0.27287, "m1": 0.34109, "m2": 0.42636, "m3": 0.53295, "m4": 0.66619, "m5": 0.83273}
    if beta == "4.47":
        #0.210476 0.263095 0.328869 0.4110859 0.5138574 0.6423218
        heavymap = {"m0": 0.210476, "m1": 0.263095, "m2": 0.328869, "m3": 0.4110859, "m4": 0.5138574, "m5": 0.6423218}
        print heavymap, heavytype
    try:
        return heavymap[heavytype]
    except:
        return None



heavy_correction = {}
heavy_correction[0.86001] = {1: 4.30584, 2: 2.60719, 3: 2.09027, 4:
                             1.83352, 5: 1.67657, 6: 1.56929, 7:
                             1.49068, 8: 1.43028, 9: 1.38226, 10:
                             1.34307, 11: 1.31044, 12: 1.2828, 13:
                             1.25909, 14: 1.23851, 15: 1.22048, 16:
                             1.20455, 17: 1.19038, 18: 1.17769, 19:
                             1.16627, 20: 1.15593, 21: 1.14654, 22:
                             1.13797, 23: 1.13013, 24: 1.12293, 25:
                             1.11629, 26: 1.11016, 27: 1.10448, 28:
                             1.09921, 29: 1.0943, 30: 1.08973, 31:
                             1.08546, 32: 1.08147, 33:1.07773,
                             34:1.07422, 35:1.07091, 36:1.06781,
                             37:1.06491, 38:1.06205, 39:1.05958,
                             40:1.05709, 41:1.05474, 42:1.05321,
                             43:1.04996, 44:1.05299, 45:1.04922,
                             46:1.05457, 47:1.06641, 48:1.06982 }



heavy_correction[0.68808] = {1: 1.98812, 2: 1.45475, 3: 1.27436, 4:
                             1.18616, 5: 1.13474, 6: 1.10159, 7:
                             1.07878, 8: 1.06236, 9: 1.05016, 10:
                             1.04086, 11: 1.03362, 12: 1.0279, 13:
                             1.02332, 14: 1.01962, 15: 1.01659, 16:
                             1.01409, 17: 1.01202, 18: 1.01029, 19:
                             1.00884, 20: 1.00762, 21: 1.00658, 22:
                             1.0057, 23: 1.00495, 24: 1.00431, 25:
                             1.00375, 26: 1.00328, 27: 1.00287, 28:
                             1.00251, 29: 1.0022, 30: 1.00193, 31:
                             1.0017, 32: 1.0015, 33:1.00132,
                             34:1.00116, 35:1.00103, 36:1.00091,
                             37:1.00081, 38:1.00071, 39:1.00063,
                             40:1.00058, 41:1.00048, 42:1.00047,
                             43:1.00037, 44:1.00055, 45:1.00046,
                             46:1.00021, 47:1.00051, 48:1.00035}

heavy_correction[0.55046] = {1: 1.51102, 2: 1.21576, 3: 1.11641, 4:
                             1.07043, 5: 1.04552, 6: 1.03071, 7:
                             1.02135, 8: 1.01518, 9: 1.01098, 10:
                             1.00806, 11: 1.00598, 12: 1.00448, 13:
                             1.00338, 14: 1.00257, 15: 1.00197, 16:
                             1.00151, 17: 1.00117, 18: 1.00091, 19:
                             1.0007, 20: 1.00055, 21: 1.00043, 22:
                             1.00034, 23: 1.00027, 24: 1.00021, 25:
                             1.00017, 26: 1.00013, 27: 1.0001, 28:
                             1.00008, 29: 1.00007, 30: 1.00005, 31:
                             1.00004, 32: 1.00003, 33:1.00003,
                             34:1.00002, 35:1.00002, 36:1.00001,
                             37:1.00001, 38:1.00001, 39:1.00001,
                             40:1.00001, 41:1., 42:1., 43:1., 44:1.,
                             45:1., 46:1., 47:1., 48:1.}

heavy_correction[0.44037] = {1: 1.33561, 2: 1.12983, 3: 1.06345, 4:
                             1.03471, 5: 1.0203, 6: 1.0124, 7:
                             1.00781, 8: 1.00504, 9: 1.00331, 10:
                             1.00221, 11: 1.00149, 12: 1.00101, 13:
                             1.0007, 14: 1.00048, 15: 1.00034, 16:
                             1.00023, 17: 1.00017, 18: 1.00012, 19:
                             1.00008, 20: 1.00006, 21: 1.00004, 22:
                             1.00003, 23: 1.00002, 24: 1.00002, 25:
                             1.00001, 26: 1.00001, 27: 1.00001, 28:
                             1., 29: 1., 30: 1., 31: 1., 32: 1., 33:
                             1, 34: 1, 35: 1, 36:1, 37:1, 38:1, 39:1,
                             40:1, 41:1, 42:1, 43:1, 44:1, 45:1,
                             46:1, 47:1, 48:1 }

# 4.35
heavy_correction[0.83273] = {1:3.59075, 2:2.25435, 3:1.83652,
                             4:1.62877, 5:1.50232, 6:1.41639, 7:1.35384, 8:1.30611, 9:1.26843,
                             10:1.2379, 11:1.21265, 12:1.19142, 13:1.17334, 14:1.15777, 15:1.14422,
                             16:1.13234, 17:1.12184, 18:1.11252, 19:1.10418, 20:1.0967, 21:1.08995,
                             22:1.08383, 23:1.07828, 24:1.07321, 25:1.06858, 26:1.06433,
                             27:1.06042, 28:1.05682, 29:1.0535, 30:1.05043, 31:1.04758, 32:1.04493,
                             33:1.04247, 34:1.04018, 35:1.03803, 36:1.03604, 37:1.03418,
                             38:1.03237, 39:1.03089, 40:1.02927, 41:1.02791, 42:1.02659,
                             43:1.02481, 44:1.02656, 45:1.02319, 46:1.02396, 47:1.03577,
                             48:1.03125}

heavy_correction[0.66619] = {1:1.8789, 2:1.39984, 3:1.23723,
                             4:1.15822, 5:1.1126, 6:1.08351, 7:1.06372, 8:1.04966, 9:1.03933,
                             10:1.03155, 11:1.02558, 12:1.02091, 13:1.01723, 14:1.01428, 15:1.0119,
                             16:1.00997, 17:1.00838, 18:1.00708, 19:1.00599, 20:1.00509,
                             21:1.00434, 22:1.00371, 23:1.00317, 24:1.00272, 25:1.00234,
                             26:1.00202, 27:1.00174, 28:1.0015, 29:1.0013, 30:1.00113, 31:1.00098,
                             32:1.00085, 33:1.00074, 34:1.00064, 35:1.00056, 36:1.00049,
                             37:1.00043, 38:1.00037, 39:1.00032, 40:1.0003, 41:1.00024, 42:1.00024,
                             43:1.0002, 44:1.00018, 45:1.00014, 46:1.00008, 47:1.00033, 48:1.00002}

heavy_correction[0.53295] = {1:1.47529, 2:1.19808, 3:1.10525,
                             4:1.0627, 5:1.0399, 6:1.02651, 7:1.01816, 8:1.01272, 9:1.00907,
                             10:1.00656, 11:1.00479, 12:1.00354, 13:1.00263, 14:1.00197,
                             15:1.00149, 16:1.00113, 17:1.00086, 18:1.00066, 19:1.0005, 20:1.00039,
                             21:1.0003, 22:1.00023, 23:1.00018, 24:1.00014, 25:1.00011, 26:1.00009,
                             27:1.00007, 28:1.00005, 29:1.00004, 30:1.00003, 31:1.00003,
                             32:1.00002, 33:1.00002, 34:1.00001, 35:1.00001, 36:1.00001,
                             37:1.00001, 38:1., 39:1., 40:1., 41:1., 42:1., 43:1., 44:1., 45:1.,
                             46:1., 47:0.999998, 48:0.999996}

heavy_correction[0.42636] = {1:1.31994, 2:1.1223, 3:1.05899,
                             4:1.03185, 5:1.01838, 6:1.01108, 7:1.00689, 8:1.00439, 9:1.00285,
                             10:1.00187, 11:1.00125, 12:1.00084, 13:1.00057, 14:1.00039,
                             15:1.00027, 16:1.00018, 17:1.00013, 18:1.00009, 19:1.00006,
                             20:1.00004, 21:1.00003, 22:1.00002, 23:1.00002, 24:1.00001,
                             25:1.00001, 26:1.00001, 27:1., 28:1., 29:1., 30:1., 31:1., 32:1.,
                             33:1., 34:1., 35:1., 36:1., 37:1., 38:1., 39:1., 40:1., 41:1., 42:1.,
                             43:1., 44:1., 45:1., 46:1., 47:1., 48:1.}

heavy_correction[0.34109] = {1:1.24529, 2:1.08691, 3:1.03865,
                             4:1.01922, 5:1.01022, 6:1.00568, 7:1.00326, 8:1.00191, 9:1.00115,
                             10:1.0007, 11:1.00043, 12:1.00027, 13:1.00017, 14:1.00011, 15:1.00007,
                             16:1.00004, 17:1.00003, 18:1.00002, 19:1.00001, 20:1.00001, 21:1.,
                             22:1., 23:1., 24:1., 25:1., 26:1., 27:1., 28:1., 29:1., 30:1., 31:1.,
                             32:1., 33:1., 34:1., 35:1., 36:1., 37:1., 38:1., 39:1., 40:1., 41:1.,
                             42:1., 43:1., 44:1., 45:1., 46:1., 47:1., 48:1.}

heavy_correction[0.27287] = {1:1.20438, 2:1.06794, 3:1.02827,
                             4:1.01315, 5:1.00654, 6:1.0034, 7:1.00182, 8:1.001, 9:1.00056,
                             10:1.00032, 11:1.00018, 12:1.00011, 13:1.00006, 14:1.00004,
                             15:1.00002, 16:1.00001, 17:1.00001, 18:1., 19:1., 20:1., 21:1., 22:1.,
                             23:1., 24:1., 25:1., 26:1., 27:1., 28:1., 29:1., 30:1., 31:1., 32:1.,
                             33:1., 34:1., 35:1., 36:1., 37:1., 38:1., 39:1., 40:1., 41:1., 42:1.,
                             43:1., 44:1., 45:1., 46:1., 47:1., 48:1.}

# 4.47
heavy_correction[0.6423218] = {1:1.77808, 2:1.3492, 3:1.20334,
                               4:1.13304, 5:1.0929, 6:1.06763, 7:1.05068, 8:1.0388, 9:1.0302,
                               10:1.02381, 11:1.01897, 12:1.01526, 13:1.01236, 14:1.01008,
                               15:1.00826, 16:1.00681, 17:1.00563, 18:1.00468, 19:1.0039, 20:1.00326,
                               21:1.00273, 22:1.0023, 23:1.00194, 24:1.00164, 25:1.00138, 26:1.00117,
                               27:1.001, 28:1.00085, 29:1.00072, 30:1.00062, 31:1.00053, 32:1.00045,
                               33:1.00039, 34:1.00033, 35:1.00028, 36:1.00024, 37:1.00021,
                               38:1.00018, 39:1.00015, 40:1.00014, 41:1.00011, 42:1.00011,
                               43:1.00009, 44:1.00008, 45:1.00009, 46:1.00005, 47:1.00008,
                               48:1.00002, 49:1., 50:1., 51:1., 52:1., 53:1., 54:1., 55:1., 56:1.,
                               57:1., 58:1., 59:1., 60:1., 61:1., 62:1., 63:1., 64:1.}

heavy_correction[0.5138574] = {1:1.44028, 2:1.18083, 3:1.09448,
                               4:1.05532, 5:1.03461, 6:1.02261, 7:1.01523, 8:1.01049, 9:1.00736,
                               10:1.00523, 11:1.00377, 12:1.00274, 13:1.002, 14:1.00148, 15:1.0011,
                               16:1.00082, 17:1.00061, 18:1.00046, 19:1.00035, 20:1.00026, 21:1.0002,
                               22:1.00015, 23:1.00012, 24:1.00009, 25:1.00007, 26:1.00005,
                               27:1.00004, 28:1.00003, 29:1.00002, 30:1.00002, 31:1.00001,
                               32:1.00001, 33:1.00001, 34:1.00001, 35:1.00001, 36:1., 37:1., 38:1.,
                               39:1., 40:1., 41:1., 42:1., 43:1., 44:1., 45:1., 46:1., 47:1., 48:1.,
                               49:1., 50:1., 51:1., 52:1., 53:1., 54:1., 55:1., 56:1., 57:1., 58:1.,
                               59:1., 60:1., 61:1., 62:1., 63:1., 64:1.}

heavy_correction[0.4110859] = {1:1.3041, 2:1.11472, 3:1.05454,
                               4:1.02902, 5:1.01651, 6:1.00981, 7:1.00602, 8:1.00378, 9:1.00241,
                               10:1.00157, 11:1.00103, 12:1.00068, 13:1.00046, 14:1.00031,
                               15:1.00021, 16:1.00014, 17:1.0001, 18:1.00007, 19:1.00005, 20:1.00003,
                               21:1.00002, 22:1.00002, 23:1.00001, 24:1.00001, 25:1.00001, 26:1.,
                               27:1., 28:1., 29:1., 30:1., 31:1., 32:1., 33:1., 34:1., 35:1., 36:1.,
                               37:1., 38:1., 39:1., 40:1., 41:1., 42:1., 43:1., 44:1., 45:1., 46:1.,
                               47:1., 48:1., 49:1., 50:1., 51:1., 52:1., 53:1.0, 54:1.0, 55:1.,
                               56:1., 57:1.0, 58:1., 59:1.0, 60:1.0, 61:1., 62:1., 63:1., 64:1.}

heavy_correction[0.328869] = {1:1.23695, 2:1.08301, 3:1.03648,
                              4:1.01793, 5:1.00942, 6:1.00517, 7:1.00293, 8:1.0017, 9:1.00101,
                              10:1.0006, 11:1.00037, 12:1.00023, 13:1.00014, 14:1.00009, 15:1.00005,
                              16:1.00003, 17:1.00002, 18:1.00001, 19:1.00001, 20:1.00001, 21:1.,
                              22:1., 23:1., 24:1., 25:1., 26:1., 27:1., 28:1., 29:1., 30:1., 31:1.,
                              32:1., 33:1., 34:1., 35:1., 36:1., 37:1., 38:1., 39:1., 40:1., 41:1.,
                              42:1., 43:1., 44:1., 45:1., 46:1., 47:1., 48:1., 49:1., 50:1., 51:1.,
                              52:1., 53:1., 54:1., 55:1., 56:1., 57:1., 58:1., 59:1., 60:1., 61:1.,
                              62:1., 63:1., 64:1.}

heavy_correction[0.263095] = {1:1.19953, 2:1.06572, 3:1.02708,
                              4:1.01247, 5:1.00614, 6:1.00316, 7:1.00168, 8:1.00091, 9:1.00051,
                              10:1.00029, 11:1.00016, 12:1.00009, 13:1.00005, 14:1.00003,
                              15:1.00002, 16:1.00001, 17:1.00001, 18:1., 19:1., 20:1., 21:1., 22:1.,
                              23:1., 24:1., 25:1., 26:1., 27:1., 28:1., 29:1., 30:1., 31:1., 32:1.,
                              33:1., 34:1., 35:1., 36:1., 37:1., 38:1., 39:1., 40:1., 41:1., 42:1.,
                              43:1., 44:1., 45:1., 46:1., 47:1., 48:1., 49:1., 50:1., 51:1., 52:1.,
                              53:1., 54:1., 55:1., 56:1., 57:1., 58:1., 59:1., 60:1., 61:1., 62:1.,
                              63:1., 64:1.}

heavy_correction[0.210476] = {1:1.17692, 2:1.0554, 3:1.02167,
                              4:1.00947, 5:1.00442, 6:1.00216, 7:1.00109, 8:1.00056, 9:1.0003,
                              10:1.00016, 11:1.00009, 12:1.00005, 13:1.00003, 14:1.00001,
                              15:1.00001, 16:1., 17:1., 18:1., 19:1., 20:1., 21:1., 22:1., 23:1.,
                              24:1., 25:1., 26:1., 27:1., 28:1., 29:1., 30:1., 31:1., 32:1., 33:1.,
                              34:1., 35:1., 36:1., 37:1., 38:1., 39:1., 40:1., 41:1., 42:1., 43:1.,
                              44:1., 45:1., 46:1., 47:1., 48:1., 49:1., 50:1., 51:1., 52:1., 53:1.,
                              54:1., 55:1., 56:1., 57:1., 58:1., 59:1., 60:1., 61:1., 62:1., 63:1.,
                              64:1.}

def get_heavy_correction(beta, heavytype):

    return heavy_correction[get_heavyq_mass(beta, heavytype)]
