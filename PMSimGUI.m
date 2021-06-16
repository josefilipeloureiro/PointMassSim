function varargout = PMSimGUI(varargin)
% PMSIMGUI MATLAB code for PMSimGUI.fig
%      PMSIMGUI, by itself, creates a new PMSIMGUI or raises the existing
%      singleton*.
%
%      H = PMSIMGUI returns the handle to a new PMSIMGUI or the handle to
%      the existing singleton*.
%
%      PMSIMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PMSIMGUI.M with the given input arguments.
%
%      PMSIMGUI('Property','Value',...) creates a new PMSIMGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PMSimGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PMSimGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PMSimGUI

% Last Modified by GUIDE v2.5 23-Apr-2018 21:58:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PMSimGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @PMSimGUI_OutputFcn, ...
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


% --- Executes just before PMSimGUI is made visible.
function PMSimGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PMSimGUI (see VARARGIN)

handles.Vel_PM = varargin{1};
handles.Deltatime_PM = varargin{2};
handles.MechanicalEnergy_PM = varargin{3};

plot((handles.Vel_PM)*3.6);
xlim([0 length(handles.Vel_PM)])
grid on
xlabel('Track Length (m)')
ylabel('Velocity (km/h)')

LapTimeText_Callback(hObject, eventdata, handles);
TopSpeedText_Callback(hObject, eventdata, handles);
EnergyConsumptionText_Callback(hObject, eventdata, handles);

% Choose default command line output for PMSimGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PMSimGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PMSimGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function LapTimeText_Callback(hObject, eventdata, handles)
% hObject    handle to LapTimeText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LapTimeText as text
%        str2double(get(hObject,'String')) returns contents of LapTimeText as a double
laptime = round(handles.Deltatime_PM(end)*10^3)/10^3;
set(handles.LapTimeText,'string',num2str(laptime));

% --- Executes during object creation, after setting all properties.
function LapTimeText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LapTimeText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function TopSpeedText_Callback(hObject, eventdata, handles)
% hObject    handle to TopSpeedText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TopSpeedText as text
%        str2double(get(hObject,'String')) returns contents of TopSpeedText as a double
topspeed = round((max(handles.Vel_PM)*3.6)*10^2)/10^2;
set(handles.TopSpeedText,'string',num2str(topspeed));

% --- Executes during object creation, after setting all properties.
function TopSpeedText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TopSpeedText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EnergyConsumptionText_Callback(hObject, eventdata, handles)
% hObject    handle to EnergyConsumptionText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EnergyConsumptionText as text
%        str2double(get(hObject,'String')) returns contents of EnergyConsumptionText as a double
energyconsumption = round((handles.MechanicalEnergy_PM(end))*10^2)/10^2;
set(handles.EnergyConsumptionText,'string',num2str(energyconsumption));

% --- Executes during object creation, after setting all properties.
function EnergyConsumptionText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EnergyConsumptionText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
