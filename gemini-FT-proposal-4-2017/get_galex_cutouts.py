import sys
sys.path.append('../../../projects')
import cutouts
import hhana
from astropy.table import Table

galex = Table.read('data/gemini-ft-sample.csv')

targets = [139, 118, 230, 27]

for ID in targets:
    prefix = 'figures/hugs-'+str(ID)
    ra, dec = galex[galex['hugs-id']==ID]['MatchRA', 'MatchDEC'][0]
    cutouts.galex.getGalexCutout(
        ra, dec, size=35, name=prefix, label_survey=False)
