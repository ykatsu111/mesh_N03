#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-


"""
Generate meshed data of the Administrative boundary from GeoJSON file.
The output file is numpy npz file.
http://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N03-v2_2.html
"""


__version__ = "1.0"


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

    if npz["lons"].ndim != 1:
        raise ValueError("lons must be 1-dim array. not %d-dim." % npz["lons"].ndim)
    if npz["lats"].ndim != 1:
        raise ValueError("lats must be 1-dim array. not %d-dim." % npz["lats"].ndim)
    if npz["data"].ndim != 2:
        raise ValueError("data must be 2-dim array. not %d-dim." % npz["data"].ndim)

    if npz["data"].shape != (npz["lats"].size, npz["lons"].size):
        raise ValueError(
            "shape of data and that of (lons, lats) are not matched.\n" +
            "shape of data: %s\n" % str(npz["data"].shape) +
            "(lons, lats) : (%d, %d)" % (npz["lats"].size, npz["lons"].size)
        )

    return npz


def get_admin_name(geoj, level):
    name = []
    for feature in geoj:
        if feature.properties[level] not in name:
            name.append(feature.properties[level])
    return name


def _get_lonlats(feature):
    import numpy as np

    if feature.geometry.type != u"Polygon":
        raise RuntimeError("Only Polygon type is supported! not %s" % feature.geometry.type)

    _coords = np.array(
        feature.geometry.coordinates
    )
    assert _coords.shape[0] == 1

    lons = _coords[0, :, 0]
    lats = _coords[0, :, 1]

    return lons, lats


def _is_surrounded(lon, lat, lons, lats):
    import numpy as np

    if lon < lons.min() or lons.max() < lon:
        return False
    if lat < lats.min() or lats.max() < lat:
        return False

    in_lon = np.array((lons <= lon) & (lon <= lons), dtype=bool)
    in_lat = np.array((lats <= lat) & (lat <= lats), dtype=bool)
    in_mesh = np.array(in_lon & in_lat, dtype=bool)

    return np.max(in_mesh) > 0


def get_surrounded_name(lon, lat, geoj, level):
    import numpy as np

    for feature in geoj:
        lons, lats = _get_lonlats(feature)
        if _is_surrounded(lon, lat, lons, lats):
            return feature.properties[level]

    raise NotFoundAdminName("Not found the surrounded administrative name!")


def _nearest_distance(lon, lat, lons, lats):
    import numpy as np

    dlons = np.abs(lons - lon)
    dlats = np.abs(lats - lat)
    dd = np.sqrt(dlons ** 2 + dlats ** 2)
    return dd.min()


def get_nearest_name(lon, lat, geoj, level):
    import numpy as np

    nearesets = np.full([len(geoj)], -np.inf)
    for i, feature in enumerate(geoj):
        lons, lats = _get_lonlats(feature)
        nearesets[i] = _nearest_distance(lon, lat, lons, lats)

    feature_min = geoj[nearesets.argmin()]
    return feature_min.properties[level]


def convert(fjson, fmask, level, msg=False, use_nearest=True):
    """
    convert GeoJSON of the Administrative boundary downloaded from nlftp.mlit.go.jp to meshed data

    :param fjson: file name of GeoJSON
    :param fmask: file name of land/sea mask data as numpy npz file
                  The file should contain following fields at least,
                  lons: 1-dimensional array presents the grid points for longitude axis
                  lats: 1-dimensional array presents the grid points for latitude axis
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
    :return: {
            'lons': 1-dimensional array presents the grid points for longitude axis,
            'lats': 1-dimensional array presents the grid points for latitude axis,
            'data': 2-dimensional array presents the administrative by the index number (0, 1, 2, ...),
            'name': the administrative name included at its corresponding index number of data
        }
    """

    import pygeoj
    import numpy as np
    from progressbar import ProgressBar
    import warnings

    # =========================================================== #
    # prepare in-data
    geoj = pygeoj.load(fjson)
    _mask = load_mask(fmask)
    mask = _mask["data"]
    lons, lats = np.meshgrid(_mask["lons"], _mask["lats"])
    # =========================================================== #

    # =========================================================== #
    # prepare out-data
    mesh = np.full_like(mask, -1, dtype=int)
    name = get_admin_name(geoj, level)
    if len(name) == 0:
        raise RuntimeError("GeoJSON has not administrative name!")
    if msg:
        print "TARGET ADMINISTRATIVE NAMES"
        for i, n in enumerate(name):
            print "\t%d: %s" % (i, n)
    # =========================================================== #

    # =========================================================== #
    # prepare progress-bar
    progress_i = 0
    progress_max = np.sum(mask)
    if msg:
        progress = ProgressBar(maxval=progress_max).start()
    # =========================================================== #

    # =========================================================== #
    # convert
    for j in range(mask.shape[0]):
        for i in range(mask.shape[1]):
            if mask[j, i]:
                if msg:
                    progress.update(progress_i + 1)

                lon = lons[j, i]
                lat = lats[j, i]

                try:
                    n = get_surrounded_name(lon, lat, geoj, level)
                except NotFoundAdminName as e:
                    if use_nearest:
                        n = get_nearest_name(lon, lat, geoj, level)
                        # warnings.warn(
                        #     "Not found the administrative name with the surrounded algorithm.\n" +
                        #     "Alternatively, the name is searched with the nearest algorithm.\n" +
                        #     "(X, Y) = (%d, %d), name=%s" % (i, j, n),
                        #     category=UserWarning
                        # )
                    else:
                        continue

                mesh[j, i] = name.index(n)

                progress_i += 1
            else:
                continue
    if msg:
        progress.finish()
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
python mesh_N03.py -j {GeoJSON} -m {Seal/Land mask} -l {AdminLevel} -o {output} [--disable-nearest]
""",
        description="""\
This program convert the administrative boundary data (GeoJSON file) 
downloaded from http://nlftp.mlit.go.jp
to meshed data.

The output file is numpy npz file containing following fields,
lons: 1-dimensional array presents the grid points for longitude axis
lats: 1-dimensional array presents the grid points for latitude axis
data: 2-dimensional array presents the administrative by the index number (0, 1, 2, ...)
name: the administrative name included at its corresponding index number of data""",
        epilog="(c) ykatsu111 (2018)",
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
lons: 1-dimensional array presents the grid points for longitude axis
lats: 1-dimensional array presents the grid points for latitude axis
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

    args = parser.parse_args()

    mesh = convert(
        args.json, args.mask, args.level,
        msg=True, use_nearest=not args.disable_nearest
    )
    np.savez(
        args.output,
        **mesh
    )
