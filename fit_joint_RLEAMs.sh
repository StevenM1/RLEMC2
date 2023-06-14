for ROI in GPe GPi PAG PPN RN SN STN Tha VTA Str
    do
    Rscript --no-save --no-restore --verbose makeSamplers_standard_singleROI.R ${ROI};
    Rscript --no-save --no-restore --verbose fit_joint_fMRI\-RLEAMs.R "./samples/dataset-trondheim_model-rleam_roi-${ROI}.RData"
done
