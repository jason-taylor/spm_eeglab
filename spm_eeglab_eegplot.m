function [fmain,fchild] = spm_eeglab_eegplot(S)
% Plot SPM MEEG data using EEGLAB's eegplot.
%  FORMAT: [fmain,fchild] = spm_eeglab_eegplot(S)
%  INPUT: Struct 'S' with fields:
%   S.D        - MEEG object or filename of MEEG object
%   S.ICA      - ICA struct (if S.compinds used) (from spm_eeglab_runica)
%   - plot options:
%   S.spacing  - Scale of main axis (def: 100)
%   S.title    - Title of main axis (def: '')
%   S.locfile  - Filename of eeglab chanloc file (default: create one)
%   - data: ONE of:
%   S.chantype - Channel type(s) to plot in main window
%   S.chaninds - Channel indices to plot in main window
%   S.chanlabs - Channel labels to plot in main window
%   S.compinds - ICA component indices (S.ICA required)
%   S.data     - data matrix (chans x timepoints [x epochs])
%   - overlay data? ONE of:
%   S.D2       - MEEG data file to superimpose
%   S.data2    - data matrix to superimpose
%   - child: if linked 'child' window desired, ONE of: 
%   S.child.chantype - as above, for child window
%   S.child.chaninds - as above, for child window
%   S.child.chanlabs - as above, for child window
%   S.child.compinds - as above, for child window
%   S.child.data     - as above, for child window
%   - child plot options:
%   S.child.spacing  - as above, for child window
%   S.child.title    - as above, for child window
%   S.child.locfile  - Filename of eeglab chanloc file (def: create one)
%
%  OUTPUT: 
%   fmain      - Figure handle of main window
%   fchild     - Figure handle of child window
%
%  spm_icameeg and spm_eeglab tools
%  by Jason Taylor (09/Mar/2018) jason.taylor@manchester.ac.uk

% Could add more options: 
% - time window
% - eegplot options

%-------------------------------------------------------------------------

fprintf('\n\n');
fprintf('++ %s\n',datestr(now));
fprintf('++ RUNNING spm_eeglab_eegplot\n')

%% Check inputs:
try data    = S.data;    catch, data    =[];   end
try data2   = S.data2;   catch, data2   =[];   end
try mytitle = S.title;   catch, mytitle = '';  end
try spacing = S.spacing; catch, spacing = 100; end
try locfile = S.locfile; catch, locfile = 0;  end

if ~isempty(data)
    ct=[]; cl=[]; ci=[]; % because we don't know channel indices, etc.
else
    try ct = S.chantype; catch, ct={}; end
    try cl = S.chanlabs; catch, cl={}; end
    try ci = S.chaninds; catch, ci=[]; end
    try ii = S.compinds; catch, ii=[]; end
end

do_chwin = isfield(S,'child');
if do_chwin
    try chdata = S.child.data; catch,chdata=[]; end
    if ~isempty(chdata)
        chct={}; chcl={}; chci=[]; chii=[];
    else
        try chct = S.child.chantype; catch, chct={}; end
        try chcl = S.child.chanlabs; catch, chcl={}; end
        try chci = S.child.chaninds; catch, chci=[]; end
        try chii = S.child.compinds; catch, chii=[]; end
    end
    try chspacing = S.child.spacing; catch, chspacing = 200; end
    try chtitle   = S.child.title;   catch, chtitle   = '';  end 
    try chlocfile = S.child.locfile; catch, chlocfile = 0;  end

end

%% Load SPM-format data file, extract data for main plot:
try D   = spm_eeg_load(S.D);  catch, error('S.D required!'); end
try D2  = spm_eeg_load(S.D2); catch, D2  = [];               end
if isfield(S,'ICA')
    if ~isstruct(S.ICA)
        load(S.ICA);
    else
        ICA = S.ICA;
        S.ICA = ICA.fname; % save RAM
    end
    timewin = ICA.timewin;
else
    ICA = [];
    timewin = [];
end

if isempty(data)
    if ~isempty(ii)
        data = ICA.activations(ii,:);
    else
        if ~isempty(ci)
            cl = chanlabels(D,ci);
        elseif ~isempty(cl)
            ci = indchannel(D,cl);
        elseif ~isempty(ct)
            cl = chanlabels(D,indchantype(D,ct));
            ci = indchannel(D,cl);
        end
        if ~iscell(timewin)
            data = selectdata(D,cl,timewin,[]);
        else
            data = [];
            for tw=1:2:size(twin)
                data = [data selectdata(D,cl,timewin{tw},[])];
            end
        end
    end
end

if ~strcmpi('continuous',D.type) && ndims(data)==2
    data = reshape(data,size(data,1),size(D,2),size(D,3));
end

%% Overlay other set of data?:
if isempty(data2)
    if ~isempty(D2)
        if ~iscell(timewin)
            data2 = selectdata(D2,cl,timewin,[]);
        else
            data2 = [];
            for tw=1:2:size(twin)
                data2 = [data2 selectdata(D2,cl,timewin{tw},[])];
            end
        end
    end
end

if ~isempty(data2)
    if ~strcmpi('continuous',D2.type) && ndims(data2)==2
        data2 = reshape(data2,size(data2,1),size(D2,2),size(D2,3));
    end
end

%% Load child-window data (if requested):
if do_chwin
    if isempty(chdata)
        if ~isempty(chii)
            chdata = ICA.activations(ii,:);
        else
            if ~isempty(chci)
                chcl = chanlabels(D,chci);
            elseif ~isempty(chcl)
                chci = indchannel(D,chcl);
            elseif ~isempty(chct)
                chcl = chanlabels(D,indchantype(D,chct));
                chci = indchannel(D,chcl);
            end
            if ~iscell(timewin)
                chdata = selectdata(D,chcl,timewin,[]);
            else
                chdata = [];
                for tw=1:size(twin)
                    chdata = [chdata selectdata(D,chcl,timewin{tw},[])];
                end
            end
        end
    end
end

%% Channel locations file (for channel labels):
if ~locfile
    % hack to use channel/IC labels rather than indices:
    if ~isempty(ci)
        locfile = 'tmp_chanlabs.xyz';
        fid=fopen(locfile,'w');
        for i=1:length(ci)
            fprintf(fid,'%d\t0\t0\t0\t%s\n',ci(i),cl{i});
        end
        fclose(fid);
    elseif ~isempty(ii)
        locfile = 'tmp_icalabs.xyz';
        fid=fopen(locfile,'w');
        for i=1:length(ii)
            fprintf(fid,'%d\t0\t0\t0\tIC%02d\n',i,ii(i));
        end
        fclose(fid);
    end
end

if ~chlocfile    
    % hack to use channel/IC labels rather than indices:
    if ~isempty(chci)
        chlocfile = 'tmp_ch_chanlabs.xyz';
        fid=fopen(chlocfile,'w');
        for i=1:length(chci)
            fprintf(fid,'%d\t0\t0\t0\t%s\n',i,chcl{i});
        end
        fclose(fid);    
    elseif ~isempty(chii)
        chlocfile = 'tmp_ch_icalabs.xyz';
        fid=fopen(chlocfile,'w');
        for i=1:length(chii)
            fprintf(fid,'%d\t0\t0\t0\tIC%02d\n',i,chii(i));
        end
        fclose(fid);
    end
end

%% Plot:
if do_chwin

    % Plot child:
    eegplot(chdata,'srate',D.fsample,'spacing',chspacing,'eloc_file',chlocfile);
    fchild = gcf;
    axes(findobj(fchild,'Tag','eegaxis'));
    if any(chtitle)
        title(chtitle);
    end
    % Re-size and re-position:
    set(fchild,'units','normalized'); 
    fchpos = get(fchild,'position');
    set(fchild,'position',[fchpos(1) .7 fchpos(3) .22]);

    % Plot main:
    if ~isempty(data2)
        eegplot(data(:,:,:),'srate',D.fsample,'children',fchild,'dispchans',min(12,size(data,1)),'data2',data2,'spacing',spacing,'eloc_file',locfile);
    else
        eegplot(data(:,:,:),'srate',D.fsample,'children',fchild,'dispchans',min(12,size(data,1)),'spacing',spacing,'eloc_file',locfile);
    end
    fmain = gcf;
    % Re-size and re-position
    set(fmain,'units','normalized'); 
    fpos = get(fmain,'position');
    set(fmain,'position',[fpos(1) fpos(2) fpos(3) .55]);

else
    
    % Plot main:
    if ~isempty(data2)
        eegplot(data(:,:,:),'srate',D.fsample,'dispchans',min(12,size(data,1)),'data2',data2,'spacing',spacing,'eloc_file',locfile);
    else
        eegplot(data(:,:,:),'srate',D.fsample,'dispchans',min(12,size(data,1)),'spacing',spacing,'eloc_file',locfile);
    end
    fmain = gcf;
    fchild = [];

end

axes(findobj(fmain,'Tag','eegaxis'));
if any(mytitle)
    title(sprintf('%s\n',mytitle));
end

%% Delete temporary chanlocs files:
% if isstr(locfile) && strcmpi('tmp',locfile(1:3))
%     delete(locfile);
% end
% if isstr(chlocfile) && strcmpi('tmp',chlocfile(1:3))
%     delete(chlocfile);
% end
