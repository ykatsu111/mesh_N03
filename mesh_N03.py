#!/usr/bin/python3.8


"""
Generate meshed data of the Administrative boundary from GeoJSON file.
The output file is numpy npz file.
http://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N03-v2_2.html
"""


__version__ = "2.3"


AdminLevels = [
    u"N03_001",
    u"N03_002",
    u"N03_003",
    u"N03_004",
    u"N03_007"
]


class NotFoundAdminName(Exception):
    pass


def load_mask(f):
    import numpy as np

    npz = np.load(f)
    if "lons" not in npz:
        raise ValueError("mask data mush contain lons field.")
    if "lats" not in npz:
        raise ValueError("mask data mush contain lats field.")
    if "data" not in npz:
        raise ValueError("mask data mush contain data field.")

    if npz["lons"].ndim != 2:
        raise ValueError("lons must be 2-dim array. not %d-dim." % npz["lons"].ndim)
    if npz["lats"].ndim != 2:
        raise ValueError("lats must be 2-dim array. not %d-dim." % npz["lats"].ndim)
    if npz["data"].ndim != 2:
        raise ValueError("data must be 2-dim array. not %d-dim." % npz["data"].ndim)

    if npz["data"].shape != npz["lats"].shape or npz["data"].shape != npz["lons"].shape:
        raise ValueError(
            "shape of data and that of (lons, lats) are not matched."
        )

    return npz


def get_admin_name_from_properties(
    properties, level, prio_ordcity=False, inc_metropolis=False, avoid_empty003=False
):
    if prio_ordcity and level != "N03_004":
        raise RuntimeError("prio_ordcity=True is only valid wiht level=N03_004")
    if inc_metropolis and level != "N03_004":
        raise RuntimeError("inc_metropolis=True is only valid with level=N03_004")
    if avoid_empty003 and level != "N03_003":
        raise RuntimeError("avoid_empty003=True is only valid with level=N03_003")
    
    name_i = properties[level]
    if prio_ordcity:
        name_i = properties["N03_003"] if "市" in properties["N03_003"] else properties["N03_004"]
    if inc_metropolis:
        name_i = properties["N03_003"] + name_i if "郡" in properties["N03_003"] else name_i
    if avoid_empty003:
        name_i = properties["N03_004"] if len(name_i) == 0 else name_i

    return name_i


def get_admin_name(geoj, level, prio_ordcity=False, inc_metropolis=False, avoid_empty003=False):
    name = []
    for feature in geoj:
        name_i = get_admin_name_from_properties(
            feature.properties, level, 
            prio_ordcity=prio_ordcity, inc_metropolis=inc_metropolis, avoid_empty003=avoid_empty003
        )
        if name_i not in name:
            name.append(name_i)
    return name


def _get_lonlats(feature):
    import numpy as np

    if feature.geometry.type != "Polygon":
        raise RuntimeError(f"Only Polygon type is supported! not {feature.geometry.type:s}")

    for _coords in feature.geometry.coordinates:
        coords = np.array(_coords, dtype=float)
        lons = coords[:, 0]
        lats = coords[:, 1]
        yield lons, lats


def crossing_number_algorithm(x, y, xx, yy):
    """
    reference: https://www.nttpc.co.jp/technology/number_algorithm.html
    """
    import numpy as np

    cross_y = ((yy[:-1] <= y) & (y < yy[1:])) | ((yy[1:] <= y) & (y < yy[:-1]))
    vt = (y - yy[:-1]) / (yy[1:] - y)
    cross_x = x < (xx[:-1] + vt * (xx[1:] - xx[:-1]))
    cn = np.sum(cross_y & cross_x)
    return cn % 2 == 1


def _is_surrounded(lon, lat, lons, lats):

    if lon < lons.min() or lons.max() < lon:
        return False
    if lat < lats.min() or lats.max() < lat:
        return False

    return crossing_number_algorithm(lon, lat, lons, lats)


def get_surrounded_name(
    lon, lat, geoj, level, 
    prio_ordcity=False, inc_metropolis=False, avoid_empty003=False
):
    for feature in geoj:
        for lons, lats in _get_lonlats(feature):
            if _is_surrounded(lon, lat, lons, lats):
                return get_admin_name_from_properties(
                    feature.properties, level, prio_ordcity=prio_ordcity,
                    inc_metropolis=inc_metropolis, avoid_empty003=avoid_empty003
                )

    raise NotFoundAdminName("Not found the surrounded administrative name!")


def calc_distance(lon1, lat1, lon2, lat2):
    import numpy as np

    x1, y1, x2, y2 = np.deg2rad(lon1), np.deg2rad(lat1), np.deg2rad(lon2), np.deg2rad(lat2)
    R = 6378.137  # km Earth radius
    dx = np.abs(np.cos(x2) * np.cos(y2) - np.cos(x1) * np.cos(y1))
    dy = np.abs(np.sin(x2) * np.cos(y2) - np.sin(x1) * np.cos(y1))
    dz = np.abs(np.sin(y2) - np.sin(y1))
    d2 = dx**2 + dy**2 + dz**2
    theta = np.arccos(1 - d2 * 0.5)
    assert np.all(~np.isnan(theta))
    return R * theta


def _nearest_distance(lon, lat, lons, lats):
    dd = calc_distance(lon, lat, lons, lats)
    return dd.min()


def get_nearest_name(
    lon, lat, geoj, level, 
    prio_ordcity=False, inc_metropolis=False, avoid_empty003=False
):
    import numpy as np

    nearesets = np.full([len(geoj)], -np.inf)
    for i, feature in enumerate(geoj):
        coods = list(_get_lonlats(feature))
        lons = np.concatenate([lons for lons, lats in coods])
        lats = np.concatenate([lats for lons, lats in coods])
        nearesets[i] = _nearest_distance(lon, lat, lons, lats)

    feature_min = geoj[nearesets.argmin()]
    return get_admin_name_from_properties(
        feature_min.properties, level, prio_ordcity=prio_ordcity,
        inc_metropolis=inc_metropolis, avoid_empty003=avoid_empty003
    )


def convert(
    fjson, fmask, level, msg=False, use_nearest=True, 
    prio_ordcity=False, inc_metropolis=False, avoid_empty003=False
):
    """
    convert GeoJSON of the Administrative boundary downloaded from nlftp.mlit.go.jp to meshed data

    :param fjson: file name of GeoJSON
    :param fmask: file name of land/sea mask data as numpy npz file
                  The file should contain following fields at least,
                  lons: 2-dimensional array presents the grid points for longitude axis
                  lats: 2-dimensional array presents the grid points for latitude axis
                  data: 2-dimensional array present land (=1) and sea (=0)
    :param level: the administrative level,
                  N03_001: prefectures
                  N03_002: branches of prefecture (only Hokkaido)
                  N03_003: metropolis and city designated by ordinance
                  N03_004: cities, towns and villages
                  N03_007: code number of administrative division
                  *** See http://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N03-v2_2.html ***
    :param msg: print message or not; default is False
                If True, your system is needed to support displaying Japanese.
    :param use_nearest: Use the nearest neighbor algorithm
                        when a grid point is not located in the administrative region
                        or not
                        This should be False if the land/sea mask data is filled only with 1.
    :param prio_ordcity: ordinance-designated city name used if the administrative range is 
                         a part of the ordinance-designated city. 
                         This option works only with 'level=N03_004'
    :param inc_metropolis: include metropolis names into the town and village names
                           This option works only with 'level=N03_004'
    :param avoid_empty003: avoid enpty name (e.g. city name)
                           Because city name without ordinance-designated city has no information in N03_003 level,
                           the output names without this option with 'level=N03_003' will be empty for the normal cities.
                           Turning to True is recommended with 'level=N03_003'.
                           THis option works only with 'level=N03_003'.
    :return: {
            'lons': 2-dimensional array presents the grid points for longitude axis,
            'lats': 2-dimensional array presents the grid points for latitude axis,
            'data': 2-dimensional array presents the administrative by the index number (0, 1, 2, ...),
            'name': the administrative name included at its corresponding index number of data
        }
    """

    import pygeoj
    import numpy as np
    from itertools import product
    from tqdm import tqdm
    import warnings

    # =========================================================== #
    # prepare in-data
    geoj = pygeoj.load(fjson)
    _mask = load_mask(fmask)
    mask, lons, lats = _mask["data"], _mask["lons"], _mask["lats"]
    # =========================================================== #

    # =========================================================== #
    # prepare out-data
    mesh = np.full_like(mask, -1, dtype=int)
    name = get_admin_name(
        geoj, level, prio_ordcity=prio_ordcity, 
        inc_metropolis=inc_metropolis, avoid_empty003=avoid_empty003
    )
    if len(name) == 0:
        raise RuntimeError("GeoJSON does not have administrative name!")
    if msg:
        print("TARGET ADMINISTRATIVE NAMES")
        for i, n in enumerate(name):
            print(f"\t{i:d}: {n:s}")
    # =========================================================== #

    # =========================================================== #
    # prepare progress-bar
    progress_max = np.sum(mask)
    if msg:
        progress = tqdm(total=progress_max)
    # =========================================================== #

    # =========================================================== #
    # convert
    for j, i in product(range(mask.shape[0]), range(mask.shape[1])):
        if mask[j, i]:
            if msg:
                progress.update(1)

            lon = lons[j, i]
            lat = lats[j, i]

            try:
                n = get_surrounded_name(
                    lon, lat, geoj, level, prio_ordcity=prio_ordcity,
                    inc_metropolis=inc_metropolis, avoid_empty003=avoid_empty003
                )
                mesh[j, i] = name.index(n)
            except NotFoundAdminName as e:
                if use_nearest:
                    n = get_nearest_name(
                        lon, lat, geoj, level, prio_ordcity=prio_ordcity,
                        inc_metropolis=inc_metropolis, avoid_empty003=avoid_empty003
                    )
                    # warnings.warn(
                    #     "Not found the administrative name with the surrounded algorithm.\n" +
                    #     "Alternatively, the name is searched with the nearest algorithm.\n" +
                    #     "(X, Y) = (%d, %d), name=%s" % (i, j, n),
                    #     category=UserWarning
                    # )
                else:
                    continue
        else:
            continue
    # =========================================================== #

    return {
        "lons": lons,
        "lats": lats,
        "data": mesh,
        "name": np.array(name)
    }


if __name__ == '__main__':
    import numpy as np
    import argparse

    parser = argparse.ArgumentParser(
        prog="mesh_N03.py",
        usage="""\
python mesh_N03.py -j {GeoJSON} -m {Seal/Land mask} -l {AdminLevel} -o {output} [--disable-nearest] [--prio_ordcity]
""",
        description="""\
This program convert the administrative boundary data (GeoJSON file) 
downloaded from http://nlftp.mlit.go.jp
to meshed data.

The output file is numpy npz file containing following fields,
lons: 2-dimensional array presents the grid points for longitude axis
lats: 2-dimensional array presents the grid points for latitude axis
data: 2-dimensional array presents the administrative by the index number (0, 1, 2, ...)
name: the administrative name included at its corresponding index number of data""",
        epilog="(c) ykatsu111 (2021)",
        add_help=True
    )
    parser.add_argument(
        "--version",
        version="%(prog)s version {__version__:s}".format(**locals()),
        action="version",
        default=False
    )
    parser.add_argument(
        "-j", "--json",
        metavar="GeoJSON",
        help="File name of GeoJSON file",
        required=True,
    )
    parser.add_argument(
        "-m", "--mask",
        metavar="Sea/Land mask",
        help="""\
File name of land/sea mask data as numpy npz file
The file should contain following fields at least,
lons: 2-dimensional array presents the grid points for longitude axis
lats: 2-dimensional array presents the grid points for latitude axis
data: 2-dimensional array present land (=1) and sea (=0)
If you do not have the mask data, an array filled 
with 1 can be assigned as the data field.
If so, note that you should not use 
the nearest neighbor algorithm (set --disable-nearest).""",
        required=True,
    )
    parser.add_argument(
        "-l", "--level",
        metavar="AdminLevel",
        help="""\
the administrative level,
N03_001: prefectures
N03_002: branches of prefecture (only Hokkaido)
N03_003: metropolis and city designated by ordinance
N03_004: cities, towns and villages
N03_007: code number of administrative division
*** See http://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N03-v2_2.html ***""",
        choices=AdminLevels
    )
    parser.add_argument(
        "-o", "--output",
        metavar="output",
        help="Output file name as numpy npz file",
        required=True
    )
    parser.add_argument(
        "--disable-nearest",
        help="""\
Use nearest neighbor algorithm or not
This algorithm will be used 
when a grid point is not located 
in the administrative region.
You should use this option if the mask data filled with 1 is assigned.""",
        required=False,
        action="store_true",
        default=False
    )
    parser.add_argument(
        "--prio-ordcity",
        help="""\
ordinance-designated city name used 
if the administrative range is 
a part of the ordinance-designated city. 
This option works only with '-l N03_004'""",
        required=False,
        action="store_true",
        default=False
    )
    parser.add_argument(
        "--inc-metropolis",
        help="""\
metropolis names will include in the town and village names.
This option works only with '-l N03_004'""",
        required=False,
        action="store_true",
        default=False
    )
    parser.add_argument(
        "--avoid-empty003",
        help="""\
avoid enpty name (e.g. city name)
Because city name without ordinance-designated city has no information in N03_003 level,
the output names without this option with '-l N03_003' will be empty for the normal cities.
Using this option is recommended with '-l N03_003'.
This option works only with '-l N03_003'.""",
        required=False,
        action="store_true",
        default=False
    )

    args = parser.parse_args()

    mesh = convert(
        args.json, args.mask, args.level,
        msg=True, use_nearest=not args.disable_nearest,
        prio_ordcity=args.prio_ordcity,
        inc_metropolis=args.inc_metropolis,
        avoid_empty003=args.avoid_empty003
    )
    np.savez(
        args.output,
        **mesh
    )