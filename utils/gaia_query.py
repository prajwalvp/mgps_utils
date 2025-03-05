import optparse
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_source"
Gaia.ROW_LIMIT = 10


def get_gaia_sources(opts):
    """
    Get all Gaia sources in field    
    """

    ref_coord = SkyCoord(ra=opts.ref_coords.split(' ')[0], dec=opts.ref_coords.split(
        ' ')[1], unit=(u.hour, u.degree), frame='icrs')
    if opts.shape == 'C':
        radius = u.Quantity(opts.search_radius, u.deg)
        query = Gaia.cone_search_async(ref_coord, radius)
        all_sources = query.get_results()
        all_sources.pprint()
    elif opts.shape == 'S':
        width = u.Quantity(opts.search_radius, u.deg)
        height = u.Quantity(opts.search_radius, u.deg)
        query = Gaia.query_object_async(
            coordinate=ref_coord, width=width, height=height)
        query.pprint()
    else:
        raise Exception(
            "Specified shape is invalid. Input C (Cone) or S (Square)")

    # Choose only specific columns of interest
if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option('--reference_coords', type=str,
                      help='Reference coordinates', dest='ref_coords')
    parser.add_option('--search_radius', type=str,
                      help='Search radius around reference coordinates', dest='search_radius', default='0.1')
    parser.add_option('--field_shape', type=str,
                      help='Choose shape of field of view around reference coords. Options: Cone (C) or Square (S)', dest='shape', default='C')
    opts, args = parser.parse_args()

    # Get Gaia sources
    get_gaia_sources(opts)
