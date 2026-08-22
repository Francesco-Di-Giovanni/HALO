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
        ("nstok", C.c_int), ("nfreq", C.c_int), ("ncount", C.c_int),
        ("intime", C.c_int), ("nscan", C.c_int), ("write_wt", C.c_int),
        ("ncard", C.c_int), ("intbase", C.c_float), ("data_format", C.c_int),
    ]


class NamesT(C.Structure):
    _fields_ = [
        ("object", C.c_char * 16), ("instrument", C.c_char * 16),
        ("cal", C.c_char * 16), ("rp_observer", C.c_char * 16),
        ("datobs", C.c_char * 12), ("datwrit", C.c_char * 12),
        ("file", C.c_char * 256), ("datsys", C.c_char * 8),
        ("version", C.c_char * 20), ("coord", C.c_char * 8),
        ("sta", C.c_char * (ANT_MAX * 8)),
        ("feed_type", C.c_char * (2 * ANT_MAX * 2)),
        ("card", C.c_char * (MAX_CARD * 80)),
        ("if_cstok", C.c_char * (4 * MAX_IF * 2)),
        ("su_name", C.c_char * (MAX_SU * 16)),
        ("su_cal", C.c_char * (MAX_SU * 4)),
        ("fg_reason", C.c_char * (MAX_FG * 24)),
        ("nx_source", C.c_char * (MAX_NX * 16)),
        ("nx_date", C.c_char * (MAX_NX * 12)),
        ("rpfitsversion", C.c_char * 20), ("bunit", C.c_char * 16),
        ("obstype", C.c_char * 16), ("errmsg", C.c_char * 80),
    ]


param_ = ParamT.in_dll(lib, "param_")
names_ = NamesT.in_dll(lib, "names_")

names_.file = FITS_PATH.encode("ascii")

MAXBUF = 500000
vis = (C.c_float * MAXBUF)()
wgt = (C.c_float * MAXBUF)()
baseline = C.c_int(0)
ut = C.c_float(0)
u = C.c_float(0)
v = C.c_float(0)
w = C.c_float(0)
flag = C.c_int(0)
bin_ = C.c_int(0)
if_no = C.c_int(0)
sourceno = C.c_int(0)


def call_rpfitsin(jstat_val):
    jstat = C.c_int(jstat_val)
    lib.rpfitsin_(
        C.byref(jstat), vis, wgt, C.byref(baseline), C.byref(ut),
        C.byref(u), C.byref(v), C.byref(w), C.byref(flag), C.byref(bin_),
        C.byref(if_no), C.byref(sourceno)
    )
    return jstat.value


# Open file
r = call_rpfitsin(-3)
print(f"open: jstat={r}")

scan_num = 0
max_scans = 400  # safety limit

while scan_num < max_scans:
    r = call_rpfitsin(-1)  # read header
    if r == 3:
        print("End of file reached while reading header.")
        break
    if r != 0:
        print(f"Unexpected jstat={r} reading header, stopping.")
        break

    scan_num += 1
    obj = names_.object.decode(errors="replace").strip()

    # Now count data records in this scan
    n_records = 0
    while True:
        r2 = call_rpfitsin(0)
        if r2 == 0:
            n_records += 1
        elif r2 == 1:
            break  # end of scan, header of next scan encountered
        elif r2 == 3:
            print(f"Scan {scan_num}: object={obj!r} nstok={param_.nstok} nfreq={param_.nfreq} records={n_records} -- EOF")
            call_rpfitsin(1)
            print("Closed file.")
            raise SystemExit
        else:
            # Intermediate table encountered (FG/SU/syscal) - library handles internally, keep reading
            pass

    print(f"Scan {scan_num}: object={obj!r} nstok={param_.nstok} nfreq={param_.nfreq} records={n_records}")

call_rpfitsin(1)
print("Closed file.")
