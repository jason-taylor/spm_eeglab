function ICA = spm_eeglab_runica(S)
% Run ICA on SPM MEEG data using EEGLAB's runica. 
%  FORMAT: ICA = spm_eeglab_runica(S)
%  INPUT: Struct 'S' with fields:
%   S.D        - MEEG object or filename of MEEG object
%   S.ncomp    - Number of components to compute
%   S.seed     - Seed for random number generator
%   S.skipbad  - skip bad channels? (1/0, def: 1) with chantype only!
%   S.prefix   - Output file prefix (default: 'ICA%d_',S.ncomp)
%   S.timewin  - [start end] time to use (s) (def: D.time([1 end]) )
%                 (Can be cell array of time windows)
%   - ONE of:
%   S.chantype - Channel type(s) to use (e.g., 'EEG')
%   S.chaninds - Channel indices to use
%   S.chanlabs - Channel labels to use
%   S.data     - data matrix (chans x timepoints [x epochs]) to use
%   -
%  OUTPUT: 
%   ICA        - Struct containing output of runica (also written to file)
%
%  spm_icameeg and spm_eeglab tools
%  by Jason Taylor (09/Mar/2018) jason.taylor@manchester.ac.uk

%-------------------------------------------------------------------------

%% Check inputs:
try ncomp   = S.ncomp;   catch, ncomp   = []; end
try seed    = S.seed;    catch, seed    = []; end
try skipbad = S.skipbad; catch, skipbad = 1;  end
try prefix  = S.prefix;  catch, prefix  = ''; end
try timewin = S.timewin; catch, timewin = []; end

try data = S.data; catch,data=[]; end
if ~isempty(data)
    ct={}; cl={}; ci=[];
else
    try ct = S.chantype; catch,ct={}; end
    try cl = S.chanlabs; catch,cl={}; end
    try ci = S.chaninds; catch,ci=[]; end
end

%% Load SPM-format data file, extract data for runica:

D = spm_eeg_load(S.D);
[~,fstem] = fileparts(D.fname);

fprintf('\n\n');
fprintf('++ %s\n',datestr(now));
fprintf('++ RUNNING spm_eeglab_runica ON %s\n',D.fname);

if isempty(timewin)
    timewin = D.time([1 end]);
end

if isempty(data)
    if isempty(cl)
        if ~isempty(ci)
            cl = chanlabels(D,ci);
        else
            fprintf('++ USING: %s channel data\n',char(ct));
            if skipbad
                fprintf('++ NOTE: Skipping any bad channels\n');
                cl = chanlabels(D,indchantype(D,ct,'GOOD'));
            else
                cl = chanlabels(D,indchantype(D,ct));
            end
        end
    end
    if ~iscell(timewin)
        fprintf('++ USING: timewindow: %g - %g\n',timewin);
        data = selectdata(D,cl,timewin,[]);
    else
        data = [];
        for tw=1:size(timewin)
            fprintf('++ USING: timewindow: %g - %g\n',timewin{tw});
            data = [data selectdata(D,cl,timewin{tw},[])];
        end
    end
end

% Reshape if epoched:
if ~strcmpi(D.type,'continuous')
    fprintf('++ USING: Epoched data\n');
    data = reshape(data,size(data,1),size(data,2)*size(data,3));
end

% Determine number of components and output filename:
if isempty(ncomp)
    ncomp = size(data,1);
end
fprintf('++ USING: Number of components: %d\n',ncomp);
if isempty(prefix)
   prefix = sprintf('ICA%d_',ncomp);
end


%% Run ICA:

if ~isempty(seed)
    fprintf('++ USING: Random seed: %g\n',seed);
    rng('default');
    rngstate = seed;
    rng(rngstate);
else
    rngstate = rng;
end
[weights,sphere,compvars,bias,signs,lrates,activations]=runica(data,'pca',S.ncomp,'extended',1);


%% Save output:
if isempty(cl)
    cl = {'custom data'};
end
ICA = [];
ICA.fname       = sprintf('%s%s',prefix,D.fname);
%--- INPUT to runica ---
ICA.D           = D.fname;
ICA.chans       = cl;
ICA.timewin     = timewin;
ICA.ncomp       = ncomp;
ICA.seed        = seed;
%--- OUTPUT of runica ---
ICA.weights     = weights; 
ICA.sphere      = sphere;
ICA.compvars    = compvars;
ICA.bias        = bias;
ICA.signs       = signs;
ICA.lrates      = lrates;
%--- Activations (memory-mapped file) ---
datfname = sprintf('%s%s.dat',prefix,fstem);
ICA.activations = file_array(datfname,size(activations),'float64');
initialise(ICA.activations);
ICA.activations(:,:) = activations;

% Save ICA:
save(ICA.fname,'ICA');

fprintf('++ Output struct ICA saved to %s\n',ICA.fname);
fprintf('++ Component activations written to %s\n (mapped to ICA.activations)\n',datfname);

return