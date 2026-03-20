
/data/nas2/software/miniconda3/bin/conda run -n graphban python /data/nas2/database/graphban/bin/predict.py \
  --trained_model /data/nas2/database/graphban/bin/trained_models/biosnap/inductive/seed12/result/best_model_epoch_45.pth \
  --save_dir graphBan_res_biosnap.csv \
  --test_path graphBan_input.csv

/data/nas2/software/miniconda3/bin/conda run -n graphban python /data/nas2/database/graphban/bin/predict.py \
  --trained_model /data/nas2/database/graphban/bin/trained_models/kiba/inductive/seed12/result/best_model_epoch_7.pth \
  --save_dir graphBan_res_kiba.csv \
  --test_path graphBan_input.csv

/data/nas2/software/miniconda3/bin/conda run -n graphban python /data/nas2/database/graphban/bin/predict.py \
  --trained_model /data/nas2/database/graphban/bin/trained_models/bindingdb/inductive/seed12/result/best_model_epoch_13.pth \
  --save_dir graphBan_res_bindingdb.csv \
  --test_path graphBan_input.csv
