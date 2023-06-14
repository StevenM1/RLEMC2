for ROI in GPe GPi PAG PPN RN SN STN Tha VTA Str
    do Rscript --no-save --no-restore --verbose fit_joint_RLEAMs_script.R ${ROI} > "./output_files/output_${ROI}.txt" 2>&1
done
