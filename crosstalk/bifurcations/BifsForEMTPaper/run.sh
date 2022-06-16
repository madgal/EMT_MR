#!/bin/bash

python emt_noEM.py &
python emt_null.py &
python mr_noWO.py &
python mr_null.py &
python psf_G_O.py &
python psf_ia_G_O.py &
