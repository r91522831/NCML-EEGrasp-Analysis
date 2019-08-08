function varargout = time_sliding(varargin)
% TIME_SLIDING MATLAB code for time_sliding.fig
%      TIME_SLIDING, by itself, creates a new TIME_SLIDING or raises the existing
%      singleton*.
%
%      H = TIME_SLIDING returns the handle to a new TIME_SLIDING or the handle to
%      the existing singleton*.
%
%      TIME_SLIDING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TIME_SLIDING.M with the given input arguments.
%
%      TIME_SLIDING('Property','Value',...) creates a new TIME_SLIDING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before time_sliding_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to time_sliding_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help time_sliding

% Last Modified by GUIDE v2.5 08-Aug-2019 16:11:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @time_sliding_OpeningFcn, ...
                   'gui_OutputFcn',  @time_sliding_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before time_sliding is made visible.
function time_sliding_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to time_sliding (see VARARGIN)

% Choose default command line output for time_sliding
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes time_sliding wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = time_sliding_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

slide_max = get(hObject,'Max');
slide_min = get(hObject,'Min');
slide_value = get(hObject,'Value');
slide_percent = (slide_value - slide_min) / (slide_max - slide_min);

time_range_id = find(handles.time_range);
time_id = floor(slide_percent * length(time_range_id)) + 1;

handles.text3.String = num2str(handles.timerstamps(time_range_id(time_id)));
time_range_id = find(handles.time_range);
handles.text4.String = num2str(handles.timerstamps(time_range_id(1)));
handles.text5.String = num2str(handles.timerstamps(time_range_id(end)));

for b = 1:handles.ncoeff
    for i = 1:handles.nfreqband
        subplot(6, handles.nfreqband, i + (b - 1) * handles.nfreqband);
        % specify ('conv', 'on') to avoid extrapolation
        topoplot(handles.freq_banded.mu_banded(i, time_range_id(time_id), b, :), handles.misc.chanlocs, 'maplimits', [handles.cmin(1, b), handles.cmax(1, b)], 'style', 'map', 'electrodes', 'off', 'conv', 'on');
        if b == 1
            title(handles.freq_banded.rg_freq_band{1, i})
        end
        if i == handles.nfreqband
            colorbar;
        end
        if i == 1
            text(-1, 0, handles.freq_banded.coeff_name{b});
        end
    end
end


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

All_dirpath = uigetdir();
handles.misc = load(fullfile(All_dirpath, 'misc.mat'));
handles.freq_banded = load(fullfile(All_dirpath, 'result_freq_banded.mat'));
handles.nfreqband = length(handles.freq_banded.rg_freq_band);
handles.timerstamps = handles.misc.tf_times{1}(:, :, 1);
handles.freqz = squeeze(handles.misc.tf_freqs{1}(:, :, 1))';
% % % freqz = freqz(2:2:68);
% % % freqz = freqz(3:end);
handles.ntime = length(handles.timerstamps);
handles.nfreq = length(handles.freqz);
handles.nelectrode = length(handles.misc.electrodes_name);
handles.coeff_name = {'\beta_0', '\beta_1', '\beta_2', '\beta_3', '\beta_4', '\beta_5'};
handles.ncoeff = length(handles.coeff_name);

handles.coeff_trace = nan(handles.ntime, handles.ncoeff);
handles.cmax = nan(1, handles.ncoeff);
handles.cmin = nan(1, handles.ncoeff);
for b = 1:handles.ncoeff
    handles.coeff_trace(:, b) = mean(reshape(permute(squeeze(handles.freq_banded.mu_banded(:, :, b, :)), [2, 1, 3]), 200, []), 2);
    tmp_mu = handles.freq_banded.mu_banded(:, :, b, :);
    handles.cmax(:, b) = max(tmp_mu(:));
    handles.cmin(:, b) = min(tmp_mu(:));
end

handles.time_range = handles.timerstamps >= -200 & handles.timerstamps < 400;
for b = 1:handles.ncoeff
    for i = 1:handles.nfreqband
        subplot(6, handles.nfreqband, i + (b - 1) * handles.nfreqband);
        % specify ('conv', 'on') to avoid extrapolation
        topoplot(handles.freq_banded.mu_banded(i, handles.time_range(1), b, :), handles.misc.chanlocs, 'maplimits', [handles.cmin(1, b), handles.cmax(1, b)], 'style', 'map', 'electrodes', 'off', 'conv', 'on');
        if b == 1
            title(handles.freq_banded.rg_freq_band{1, i})
        end
        if i == handles.nfreqband
            colorbar;
        end
        if i == 1
            text(-1, 0, handles.freq_banded.coeff_name{b});
        end
    end
end

guidata(hObject, handles);
