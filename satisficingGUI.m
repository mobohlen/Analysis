function varargout = satisficingGUI(varargin)
% SATISFICINGGUI MATLAB code for satisficingGUI.fig
%      SATISFICINGGUI, by itself, creates a new SATISFICINGGUI or raises the existing
%      singleton*.
%
%      H = SATISFICINGGUI returns the handle to a new SATISFICINGGUI or the handle to
%      the existing singleton*.
%
%      SATISFICINGGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SATISFICINGGUI.M with the given input arguments.
%
%      SATISFICINGGUI('Property','Value',...) creates a new SATISFICINGGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before satisficingGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to satisficingGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help satisficingGUI

% Last Modified by GUIDE v2.5 19-Apr-2016 15:32:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @satisficingGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @satisficingGUI_OutputFcn, ...
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

% --- Executes just before satisficingGUI is made visible.
function satisficingGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn. 
% First function to run once open; initializes variables with defaults.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to satisficingGUI (see VARARGIN)

%% Initialize lists of files with files from "Arya" folder
% popupmenu.String: Folder name containing *.mat data files
% popupmenu.UserData: Full path name for corresponding folder
% All available folders are initialized in the String field.
% UserData only contains the path of the active folder.
dataFolders = {'Arya','Lana'};
filePath = ['../Data/',dataFolders{1}];
fileDir = dir(fullfile(filePath,'*.mat'));
fileList = cell(length(fileDir),1);
set(handles.popupmenu,'String',dataFolders,'Value',1,'UserData',filePath);
for i = 1:length(fileDir)
    fileList{i} = fileDir(i).name;
end
%% Sort files in list in descending chronological order and add to listbox
[files,allDates] = readAndSortDates(filePath,fileList);
set(handles.listLeft,'String',allDates,'Value',1,...
        'UserData',files);

%% Default Run options
set(handles.plotIndiv,'Value',1);
set(handles.plotAll,'Value',1);
set(handles.correctFlag,'Enable','on');
set(handles.xdateFlag,'enable','on');
set(handles.dispTTFlag,'enable','on');

%% Choose default command line output for satisficingGUI
handles.output = hObject;
%% Update handles structure
guidata(hObject, handles);
% UIWAIT makes satisficingGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = satisficingGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% --------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uiFile_ClickedCallback(hObject, eventdata, handles)
dir = uigetdir();
set(hObject,'UserData',dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function listLeft_Callback(hObject, eventdata, handles)
% --- Executes on selection change in listLeft. (Unused)

function listLeft_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
%%  Allow for multiple selections (up to 1000)
set(hObject,'Max',1000,'Min',0);
%% Set color
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function listRight_Callback(hObject, eventdata, handles)
% --- Executes on selection change in listRight. (Unused)

function listRight_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
%% Initialize list and allow for multiple selections (up to 1000)
set(hObject,'String',[]);
set(hObject,'Max',1000,'Min',0);
%% Set color
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushRight_Callback(hObject, eventdata, handles)
% --- Executes on button press in pushRight.
% Adds selected file from listLeft to right list for analysis.
% Adds name to "String" for display on listRight
% Adds filename to "UserData" on listRight
% Sorts UserData and String in descending chronological order
%% Get selected file on left list
selectedVals = get(handles.listLeft,'Value');
listLeft = get(handles.listLeft,'String');
fileList = get(handles.listLeft,'UserData');
strToAdd = listLeft(selectedVals);
filesToAdd = fileList(selectedVals);
%% Add selected file to right list if not already there
strRight = get(handles.listRight,'String');
filesRight = get(handles.listRight,'UserData');
for i = 1:length(strToAdd)
    % Loop through selected files
    if(~any(ismember(strRight,strToAdd{i})))
        % Add files if not already there
        rLen = length(strRight);
        strRight{rLen+1,1} = strToAdd{i};
        filesRight{rLen+1,1} = filesToAdd{i};
    end
end
%% Resort and display
[~,sInd] = sort(datenum(strRight));
set(handles.listRight,'String',flipud(strRight(sInd)),...
    'UserData',flipud(filesRight(sInd)));
%% Update selected value on listLeft
if(selectedVals<length(listLeft) & length(strToAdd)==1)
    set(handles.listLeft,'Value',selectedVals+1);
end
%%
function pushLeft_Callback(hObject, eventdata, handles)
% --- Executes on button press in pushLeft.
% Removes selected file from listRight
% Removes name from "String"  on listRight
% Removes filename from "UserData" on listRight
% Sorts UserData and String in descending chronological order
%% Get selected file on right list
selectedVals = get(handles.listRight,'Value');
strRight = get(handles.listRight,'String');
filesRight = get(handles.listRight,'UserData');
%% Remove file from right list
if(length(strRight)>0)
    strToRemove = strRight(selectedVals);
    for i = 1:length(strToRemove)
        % Loop through selected strings/files
        for j = 1:length(strRight)
            % Loop through all strings/files
            if(strcmp(strToRemove(i),strRight{j}))
                % Remove selected strings/files
                strRight(j) = [];
                filesRight(j) = [];
                break
            end
        end
    end
    % Update selected file on listRight
    if(length(strToRemove)>1)
        set(handles.listRight,'Value',1);
    else
        if(selectedVals > 1)
         set(handles.listRight,'Value',selectedVals-1)
        end
    end
end
%% Update listRight
set(handles.listRight,'String',strRight,...
    'UserData',filesRight);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function selAllLeft_Callback(hObject, eventdata, handles)
% --- Executes on button press in selAllLeft.
% Selects all elements in listLeft
list = get(handles.listLeft,'String');
set(handles.listLeft,'Value',1:length(list));

function selAllRight_Callback(hObject, eventdata, handles)
% --- Executes on button press in selAllRight.
% Selects all elements in listRight
list = get(handles.listRight,'String');
set(handles.listRight,'Value',1:length(list));

function deselectLeft_Callback(hObject, eventdata, handles)
% --- Executes on button press in deselectLeft.
% Deselects all elements in left listbox
set(handles.listLeft,'Value',[]);

function deselectRight_Callback(hObject, eventdata, handles)
% --- Executes on button press in deselectRight.
% Deselects all elements in right listbox
set(handles.listRight,'Value',[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function popupmenu_Callback(hObject, eventdata, handles)
% --- Executes on selection change in popupmenu.
% String: Folder name containing *.mat data files
% UserData: Full path name for corresponding folder
% All available folders are initialized in "satisficingGUI_OpeningFcn" in
% the String field. UserData only contains the path of the active folder.

%% Obtain list of files from selected folder
dataFolders = get(hObject,'String');
val = get(hObject,'Value');
filePath = ['../Data/',dataFolders{val}];
fileDir = dir(fullfile(filePath,'*.mat'));
fileList = cell(length(fileDir),1);
set(hObject,'UserData',filePath);
for i = 1:length(fileDir)
    fileList{i} = fileDir(i).name;
end
%% Sort files in descending chronological order and add to listLeft
[files,allDates] = readAndSortDates(filePath,fileList);
set(handles.listLeft,'String',allDates,'Value',1,...
        'UserData',files);
%% Reset listRight
set(handles.listRight,'String',[],'UserData',[]);

function popupmenu_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Radio button callbacks (currently unused)
function choiceFlag_Callback(hObject, eventdata, handles)
% --- Executes on button press in choiceFlag.
function logRegFlag_Callback(hObject, eventdata, handles)
% --- Executes on button press in logRegFlag.
function subWeightFlag_Callback(hObject, eventdata, handles)
% --- Executes on button press in subWeightFlag.
function modelFlag_Callback(hObject, eventdata, handles)
% --- Executes on button press in modelFlag.
function responseTimeFlag_Callback(hObject, eventdata, handles)
% --- Executes on button press in responseTimeFlag.
function modelCorrFlag_Callback(hObject, eventdata, handles)
% --- Executes on button press in modelCorrFlag.
function correctFlag_Callback(hObject, eventdata, handles)
% --- Executes on button press in correctFlag.
function plotIndiv_Callback(hObject, eventdata, handles)
% --- Executes on button press in plotIndiv.
function ideal_Callback(hObject, eventdata, handles)
% --- Executes on button press in ideal.

function dispTTFlag_Callback(hObject, eventdata, handles)
% --- Executes on button press in dispTTFlag.
function xdateFlag_Callback(hObject, eventdata, handles)
% --- Executes on button press in xdateFlag.


function plotAll_Callback(hObject, eventdata, handles)
% --- Executes on button press in plotAll.
selfState = get(hObject,'Value');
if(selfState == 1)
    %set(handles.correctFlag,'enable','on');
    set(handles.xdateFlag,'enable','on');
    set(handles.dispTTFlag,'enable','on');
else
    %set(handles.correctFlag,'enable','off');
    set(handles.xdateFlag,'enable','off');
    set(handles.dispTTFlag,'enable','off');
end

function closeFigs_Callback(hObject, eventdata, handles)
% --- Executes on button press in closeFigs.
% Runs function closeFigs();
closeFigs();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions

function closeFigs()
% Closes all figures besides GUI
for i = 1:1000
    try
        close(i)
    catch
    end
end

function [files,dates] = readAndSortDates(filePath,fileList)
% Perform initial reading sorting and of dates when opening new folder
% Sorts files in descending chronological order
d = 0;
for i =  1:length(fileList)
    tempFile = fileList{i};
    tempRead = load(fullfile(filePath,tempFile));
    d = d+1;
    files{d} = fileList{i};
    dateSort{d} = tempRead.header.date; %location of date in file name
end
[~,sortIndex] = sort(datenum(dateSort)); 
files = flipud(files(sortIndex)'); % reorder summary files by date
dates = flipud(dateSort(sortIndex)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function runAnalysis_Callback(hObject, eventdata, handles)
% --- Executes on button press in runAnalysis.
% Main analyses (all the real work is done here)

% Requires following functions to run:
%   cueChoice.m
%   logRegPerformance.m
%   logdet.m
%   bayes_logit_fit.m

closeFigs();
if(length(get(handles.listRight,'UserData'))>0)
    flags = [];
    % IMPORTANT: "flags" must match up with appropriate flags in
    % runAnalysis.m
    flags(1) = get(handles.plotIndiv,'Value'); % plot individual days
    flags(2) = get(handles.plotAll,'Value'); % plot cumulative data
    flags(3) = get(handles.choiceFlag,'Value'); % choice probabilities
    flags(4) = get(handles.logRegFlag,'Value'); % logisitc regression curve
    flags(5) = get(handles.subWeightFlag,'Value'); % subjective weights
    flags(6) = get(handles.modelFlag,'Value'); % Bayes model comparison
    flags(7) = get(handles.responseTimeFlag,'Value'); % response times histogram
    flags(8) = get(handles.correctFlag,'Value'); % plot proportion correct
    flags(9) = get(handles.modelCorrFlag,'Value'); % model correlation
    flags(10) = get(handles.ideal,'Value'); % indicate ideal observer performance on figures
    flags(11) = get(handles.xdateFlag,'Value'); % xlabel with date vs. number
    flags(12) = get(handles.dispTTFlag,'Value'); % display trial times on plots
    
    files = flipud(get(handles.listRight,'UserData'));
    % Note: needs to be flipped because flipped in display
    filePath = get(handles.popupmenu,'UserData');
    
    % Run analysis function for selected files,flags
    run_fullAnalysis(filePath,files,flags);

end
