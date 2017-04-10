import os
import numpy as np
import hhana
from toolbox.astro import angsep
from IPython.display import Image, display, HTML


__all__ = ['clusters', 'get_sample', 'get_nsa', 'get_yang', 
           'find_nearest_in_Mpc', 'display_candy', 'display_group']


clusters = {'3.1': {'coord': [181, 1.5], 'size': 1.0}, 
            '3.2': {'coord': [182.5, 1.], 'size': 0.5},
            '1.1': {'coord': [37.0, -2.], 'size': 1.7},
            '2.1': {'coord': [140.0, 1.], 'size': 1.2},
            '3.3': {'coord': [180.5, -0.5], 'size': 1.},
            '4.1': {'coord': [225, 1], 'size': 1.8},
            '4.2': {'coord': [211, -1], 'size': 1.5}}


def get_sample(min_size=3.2):
    sample = hhana.SuperCat('hsc-hugs')
    sample.make_cuts((sample['candy']==1) &
                     (sample['FLUX_RADIUS(i)']>min_size/0.168))
    return sample


def get_nsa(max_mass=5e10, min_z=0.001):
    nsa = hhana.SuperCat('nsa')
    nsa.make_cuts((nsa['mass_h07']>max_mass) & (nsa['z']>min_z))
    return nsa


def get_yang(logMh_min=12.7, zmin=0.01, zmax=0.05):
    yang = hhana.SuperCat('yang_groups')
    yang.make_cuts(
	(yang['logMh(Lest)'] > logMh_min) &
	(yang['z'] > zmin) &
	(yang['z'] < zmax)
    )
    return yang


def find_nearest_in_Mpc(sam, ref, cosmo=None):
    if cosmo is None:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(70, 0.3)
    ang_diam_dists = cosmo.angular_diameter_distance(ref['z'])
    ref_idx = []
    sep_Mpc = []
    for ra, dec in sam['ra', 'dec']:
        seps = angsep(ref['ra'], ref['dec'], ra, dec, sepunits='radian')
        seps *= ang_diam_dists.value
        sep_Mpc.append(seps.min())
        ref_idx.append(seps.argmin())
    ref_idx = np.array(ref_idx)
    sep_Mpc = np.array(sep_Mpc)
    return ref_idx, sep_Mpc


def display_candy(num, size=300):
    image_dir = '../udg-zoo/io/images/'
    display(
        Image(filename=image_dir+'candy-'+str(num)+'.png', 
              height=size, width=size)
    )

def display_group(ids, size=220):
    image_dir = '../udg-zoo/io/images/'
    files = ['candy-'+str(num)+'.png' for num in ids]
    for fn in files:
        new_fn = os.path.join('images', fn)
        if not os.path.isfile(new_fn):
            shutil.copyfile(os.path.join(image_dir, fn), new_fn)
    files = ['images/candy-'+str(num)+'.png' for num in ids]
    imglist=''.join(["<img style='width: "+str(size)+\
                     "px; height:"+str(size)+"px;"
                     "margin:0px; float: left; border: 5px solid black;'" 
                     "src='%s' />" % str(s) for s in files])
    display(HTML(imglist))
