import ctypes as C

LIB_PATH = "/home/franz/miniconda3/envs/rpfits_env/lib/librpfits.so"
FITS_PATH = "/home/franz/FRO/Parsers/ATCA/Surveys/2002-02-09_2302.C1025"

ANT_MAX = 16
MAX_IF = 48
MAX_SU = 2048
MAX_CARD = 2240
MAX_FG = 32
MAX_NX = 256

lib = C.CDLL(LIB_PATH)


class ParamT(C.Structure):
    _fields_ = [
        ("nstok", C.c_int),
        ("nfreq", C.c_int),
        ("ncount", C.c_int),
        ("intime", C.c_int),
        ("nscan", C.c_int),
        ("write_wt", C.c_int),
        ("ncard", C.c_int),
        ("intbase", C.c_float),
        ("data_format", C.c_int),
    ]


class NamesT(C.Structure):
    _fields_ = [
        ("object", C.c_char * 16),
        ("instrument", C.c_char * 16),
        ("cal", C.c_char * 16),
        ("rp_observer", C.c_char * 16),
        ("datobs", C.c_char * 12),
        ("datwrit", C.c_char * 12),
        ("file", C.c_char * 256),
        ("datsys", C.c_char * 8),
        ("version", C.c_char * 20),
        ("coord", C.c_char * 8),
        ("sta", C.c_char * (ANT_MAX * 8)),
        ("feed_type", C.c_char * (2 * ANT_MAX * 2)),
        ("card", C.c_char * (MAX_CARD * 80)),
        ("if_cstok", C.c_char * (4 * MAX_IF * 2)),
        ("su_name", C.c_char * (MAX_SU * 16)),
        ("su_cal", C.c_char * (MAX_SU * 4)),
        ("fg_reason", C.c_char * (MAX_FG * 24)),
        ("nx_source", C.c_char * (MAX_NX * 16)),
        ("nx_date", C.c_char * (MAX_NX * 12)),
        ("rpfitsversion", C.c_char * 20),
        ("bunit", C.c_char * 16),
        ("obstype", C.c_char * 16),
        ("errmsg", C.c_char * 80),
    ]


param_ = ParamT.in_dll(lib, "param_")
names_ = NamesT.in_dll(lib, "names_")

names_.file = FITS_PATH.encode("ascii")

jstat = C.c_int(-3)  # -3 = open file for reading
vis = (C.c_float * 4)()
wgt = (C.c_float * 4)()
baseline = C.c_int(0)
ut = C.c_float(0)
u = C.c_float(0)
v = C.c_float(0)
w = C.c_float(0)
flag = C.c_int(0)
bin_ = C.c_int(0)
if_no = C.c_int(0)
sourceno = C.c_int(0)

lib.rpfitsin_(
    C.byref(jstat), vis, wgt, C.byref(baseline), C.byref(ut),
    C.byref(u), C.byref(v), C.byref(w), C.byref(flag), C.byref(bin_),
    C.byref(if_no), C.byref(sourceno)
)

print("jstat after open call:", jstat.value)
print("object:", names_.object)
print("instrument:", names_.instrument)
print("rpfitsversion:", names_.rpfitsversion)
print("errmsg:", names_.errmsg)

# --- Read scan header (jstat = -1) ---
jstat2 = C.c_int(-1)
lib.rpfitsin_(
    C.byref(jstat2), vis, wgt, C.byref(baseline), C.byref(ut),
    C.byref(u), C.byref(v), C.byref(w), C.byref(flag), C.byref(bin_),
    C.byref(if_no), C.byref(sourceno)
)
print("---")
print("jstat after header read:", jstat2.value)
print("object:", names_.object)
print("instrument:", names_.instrument)
print("datobs:", names_.datobs)
print("rp_observer:", names_.rp_observer)
print("nstok:", param_.nstok)
print("nfreq:", param_.nfreq)
print("nscan:", param_.nscan)
print("errmsg:", names_.errmsg)

# --- Close file (jstat = 1) ---
jstat3 = C.c_int(1)
lib.rpfitsin_(
    C.byref(jstat3), vis, wgt, C.byref(baseline), C.byref(ut),
    C.byref(u), C.byref(v), C.byref(w), C.byref(flag), C.byref(bin_),
    C.byref(if_no), C.byref(sourceno)
)
print("jstat after close:", jstat3.value)
