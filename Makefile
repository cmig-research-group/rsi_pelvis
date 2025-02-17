MATLAB = /usr/pubsw/packages/matlab/R2023a
MCC = $(MATLAB)/bin/mcc
MATFLAGS = -I ../cmig_utils/B0_distortion -I ../cmig_utils/ctxmgh -I ../cmig_utils/dicom -I ../cmig_utils/dmri -I ../cmig_utils/dmri/model_fitting \
 	-I ../cmig_utils/eddy_current -I ../cmig_utils/gradient_nonlinearity -I ../cmig_utils/matlab -I ../cmig_utils/mmps \
	-I ../cmig_utils/morphometry -I ../cmig_utils/morphometry/volgetvxlsvalMEX -I ../cmig_utils/morphometry/volgetCOMMEX -I ../cmig_utils/morphometry/resliceMEX \
	-I ../cmig_utils/motion -I ../cmig_utils/dicom/dicts -I ../cmig_utils/ctxmgh/coords/ \
	-I ../cmig_utils/B0_distortion/mrpg -I utils -I preprocessing -I autoseg 

ALL_PROGS = \
	RSI_pipeline_multisort/RSI_pipeline_multisort \
	contour_prostate/contour_prostate

all: $(ALL_PROGS)

RSI_pipeline_multisort/RSI_pipeline_multisort: RSI_pipeline_multisort.m RSI_pipeline.m Makefile
	$(MCC) -m $< -d $(basename $<) $(MATFLAGS)

contour_prostate/contour_prostate: autoseg/contour_prostate.m
	$(MCC) -m $< -d $(basename $<) $(MATFLAGS)

clean:
	rm -rf RSI_pipeline_multisort

MMPS_DIR=/usr/pubsw/packages/MMPS/MMPS_267
CWD := $(shell pwd)
TEST_OUTPUT := $(CWD)/data/test_output

run: RSI_pipeline_multisort/RSI_pipeline_multisort
	cd data/ && MMPS_DIR=$(MMPS_DIR) ../RSI_pipeline_multisort/run_RSI_pipeline_multisort.sh $(MATLAB) example_images $(TEST_OUTPUT) ../example_params.m

run_docker_test: contour_prostate/contour_prostate
	cd data/ && MMPS_DIR=$(MMPS_DIR) ../autoseg/contour_prostate/run_contour_prostate.sh $(MATLAB) $(TEST_OUTPUT)/processed_outputs/PDS3412/20240105/Series_5__RSI_TE_short/T2_corrected_GUW.mgz Docker

install: autoseg/call_docker.sh ../cmig_utils/dicom/dicts/gems-dicom-dict.txt
	mkdir -p data/example_images
	cp autoseg/call_docker.sh data/
	cp ../cmig_utils/dicom/dicts/gems-dicom-dict.txt data/