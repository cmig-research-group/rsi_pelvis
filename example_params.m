% Example configuration file for processing a prostate RSI exam
% The path to this file (or similar) should be the 3rd argument to RSI_pipeline_multisort or RSI_pipeline
% Text-string parameters are case insensitive (e.g., Docker = docker)  

% RSI model parameters
params.ModelADCs = [1.0528e-04 0.0018 0.0036 0.1220];

% Explicitly declare tensor file (if not declared in DICOM header)
% params.TensorFile = 'tensor97.dat';

% B0 distortion correction 
params.B0DISCO = 1;
  % RPG options (if using RPG)
  % params.RPG.Installation = 'docker'; % Enable if RPG is installed via Docker
  params.RPG.B0optFlag = 1;
  params.RPG.B0corrmin = 0.95;

% Other corrections
params.GradWarpFlag = 1;
params.EddyCorrFlag = 1;
params.MotionCorrFlag = 1;
params.CorrectNoise = 1;

% Automated prostate segmentation
  % Supports deep-learning AI segmentation models from two vendors: CMIG (default) and Cortechs.ai
  % Both models can be run using Docker, Podman (with alias docker=podman), or Singularity (default for CMIG users)
    % If on the CMIG network, set the container to 'Singularity'
    % If you have and want to use the Docker/Podman container, set the container to 'Docker' 
params.ProstateSeg = 1;
params.ProstateSegVendor = 'CMIG'; % Use 'Cortechs' for the model from Cortechs.ai (if you have the Docker image)
params.ProstateSegContainer = 'Singularity';

% Automated urethra segmentation
  % Deep-learning AI segmentation model from CMIG
  % Can be run using Docker, Podman (with alias docker=podman), or Singularity (default for CMIG users)
    % If on the CMIG network, set the container to 'Singularity'
    % If you have and want to use the Docker/Podman container, set the container to 'Docker' 
params.UrethraSeg = 1;
params.UrethraSegContainer = 'Docker';

% DICOM outputs
params.WriteDICOMS = 1;
  params.SelectDICOMS.T2 = 1;
  params.SelectDICOMS.DCE = 1;
  params.SelectDICOMS.DWI = 0;
  params.SelectDICOMS.ADC = 1;
  params.SelectDICOMS.RSICmaps = 0;
  params.SelectDICOMS.RSIrs = 1;
  params.SelectDICOMS.RSIrsReportOverlay = 1;
  params.SelectDICOMS.ProstateAutoSeg_RT = 1;
  params.SelectDICOMS.ProstateAutoSeg_SEG = 1;
  params.SelectDICOMS.UrethraAutoSeg_RT = 1;
  params.SelectDICOMS.LesionAutoSeg_RT = 1;
  params.SelectDICOMS.LesionAutoSeg_SEG = 1;

% Python vitrual environment (only necessary to generate DICOM SEG contours)
params.PythonVEnv = '/usr';

% QC
params.DebugFlag = 0;
params.CheckProtocolCompliance = 1;
params.ProtocolReference = 'artpro_protocol_reference.mat';

% Save non-RSI series as mgz files 
params.WriteOptionalSeries = 0;
