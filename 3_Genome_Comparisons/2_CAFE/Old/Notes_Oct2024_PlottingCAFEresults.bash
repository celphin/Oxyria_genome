#####################
# Plotting CAFE results
# Oct 2024
# https://github.com/moshi4/CafePlotter
##########################

cd /lustre04/scratch/celphin/Oxyria/CAFE/

module load StdEnv/2023 python/3.11.5 scipy-stack/2023b
pip install cafeplotter

cafeplotter -i /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/output/error_model \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Cafe_plots_Brassicaceae 

cafeplotter -i /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/error_model \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Cafe_plots_Polygonaceae 

cafeplotter -i /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/output/error_model \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Cafe_plots_Rosaceae 

#-------------------------------------------
cafeplotter -i /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/output/error_model \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Cafe_plots_Brassicaceae --format pdf

cafeplotter -i /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/error_model \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Cafe_plots_Polygonaceae --format pdf

cafeplotter -i /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/output/error_model \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Cafe_plots_Rosaceae --format pdf


#########################################
# Plotting CAFE results - does not work yet
# https://github.com/LKremer/CAFE_fig
# https://github.com/etetoolkit/ete/issues/354

cd /lustre04/scratch/celphin/Oxyria/CAFE/

module load StdEnv/2023 python/3.11.5 scipy-stack/2023b qt/6.5.3 
pip3 install 'ete3==3.0.0b35'

wget https://raw.githubusercontent.com/LKremer/CAFE_fig/refs/heads/master/CAFE_fig.py

python3 /lustre04/scratch/celphin/Oxyria/CAFE/CAFE_fig.py \
/lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/error_model/Base_report.cafe \
-pb 0.05 -pf 0.05 --dump test/ -g .pdf --count_all_expansions
