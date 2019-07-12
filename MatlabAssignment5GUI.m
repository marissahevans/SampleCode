%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% MACHINE LEARNING GUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Date Created: July 7th 2019
% Date Last Edit: July 10th 2019
% Author: Marissa Evans, mhe229@nyu.edu

% Assumptions: The file 'nycMiddleSchools.xlsx' is in the same 
% directory as the code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a GUI that allows for the user the interact with the data in
% various ways and recieve different outputs depending on their choices. 

% Load Data: this button loads the data set mentioned above (if the user
% needs to start over at any point all they must do is reload the data)

% Remove Missing Data By School: This removes the NaN's row wise,
% illiminating all data for a school if some of the information is missing.

% Remove Missing Data Using Average: This uses inputation to remove the
% NaN's by replacing them with the mean for the column they are in, keeping
% all schools present. 

% Slider: This slider allows the user to choose a cut off point for schools
% that represent a given percentage of the total admissions. (i.e. 1%, 2%
% 3% etc.)

% Recode: This button changes the schools admissions percentage to either a
% 0 or a 1, depending on if they are above or below the cut off selected by
% the user. After recoding it shows how many schools are included in each
% group. 

% Run a PCA: This button runs a PCA on the data after it was recoded into
% catagorical outputs and shows the scree plot on the central axis. 

% Factor Comparison: Choose one of the 3x top factors to compare to the
% other on separate axes using the radio buttons. Click 'show figure' to
% display the figure on the central axis. 

% Cluster Analysis: If the toggle box is selected the cluster will use the
% 'optimal' choice based on a sillouhette, if it is not selected the user
% is able to input the number of factors they would like to see into the
% box. This uses the two factors selected above as the loadings. 

% Clear Figure: Between presenting different displays on the central axis,
% select the clear figure button for optimal performance and display
% accuracy. 

% Run a Classification: Select either to train/test with an 80/20 split or
% an all but 1 split. Then select the type of model you would like to use
% below (either SVM or Random Forest- if using the RF model input the
% amount of trees in the box below). The model accuracy will be output
% depending on the choices made above. 


% OUTPUTS:
%GUI opens various graphs depending on input choices, additionally outputs
%the accuracy of the model depending on the parameters selected. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = MatlabAssignment5GUI(varargin)
% MATLABASSIGNMENT5GUI MATLAB code for MatlabAssignment5GUI.fig
%      MATLABASSIGNMENT5GUI, by itself, creates a new MATLABASSIGNMENT5GUI or raises the existing
%      singleton*.
%
%      H = MATLABASSIGNMENT5GUI returns the handle to a new MATLABASSIGNMENT5GUI or the handle to
%      the existing singleton*.
%
%      MATLABASSIGNMENT5GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MATLABASSIGNMENT5GUI.M with the given input arguments.
%
%      MATLABASSIGNMENT5GUI('Property','Value',...) creates a new MATLABASSIGNMENT5GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MatlabAssignment5GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MatlabAssignment5GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MatlabAssignment5GUI

% Last Modified by GUIDE v2.5 08-Jul-2019 15:54:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @MatlabAssignment5GUI_OpeningFcn, ...
    'gui_OutputFcn',  @MatlabAssignment5GUI_OutputFcn, ...
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


% --- Executes just before MatlabAssignment5GUI is made visible.
function MatlabAssignment5GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MatlabAssignment5GUI (see VARARGIN)

% Choose default command line output for MatlabAssignment5GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MatlabAssignment5GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MatlabAssignment5GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadData.
function loadData_Callback(hObject, eventdata, handles)
% hObject    handle to loadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[data textData rawData] = xlsread('nycMiddleSchools.xlsx'); %input data from Excell format
handles.data = data; %set up a handle for the data
set(handles.cleanRowWise, 'Visible', 'On'); %turn the clean data row-wise button on
set(handles.cleanAverage, 'Visible', 'On'); %turn the clean data using imputation button on
guidata(hObject, handles);


% --- Executes on button press in cleanRowWise.
function cleanRowWise_Callback(hObject, eventdata, handles)
% hObject    handle to cleanRowWise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nanCount = sum(isnan(handles.data),2); %counts missing values per row
nanCount = find(nanCount==0) ; %finds all rows without NaNs
dataRowWise = handles.data(nanCount,:) ; %selects only the rows without NaN's in the main data stucture
handles.data = dataRowWise; %updates the data handle with new data

studentAdmin = round(handles.data(:,1)); %calculates the number of students admitted form each scool
totalAdmin = sum(handles.data(:,1)); %calculates the total number of students admitted from all schools
for ii = 1 : length(handles.data) %updates each school's number of students to the proportion of total students that came from that school
    handles.data(ii,1) = studentAdmin(ii)/totalAdmin ;
end

set(handles.accStudSlide, 'Visible', 'On'); %turn slider on
set(handles.zeroPercent,'Visible','On'); %turn left slider label on
set(handles.twoPercent,'Visible','On'); %turn center slider label on
set(handles.fourPercent,'Visible','On'); %turn right slider label on
set(handles.sliderText,'Visible','On'); %turn slider instructions on
set(handles.recode, 'Visible','On'); %turn recode button on
set(handles.sliderVal, 'Visible','On'); %turn the slider output box on
set(handles.text16, 'Visible', 'On'); %turn reload instructions on
guidata(hObject, handles);


% --- Executes on button press in cleanAverage.
function cleanAverage_Callback(hObject, eventdata, handles)
% hObject    handle to cleanAverage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cleanData = ones(length(handles.data),15)*99; %preallocate a matrix for the clean data
for ii = 1:16 %run though loop 16x for all columns
    data2 = handles.data; %make a copy of data
    columnAverage = nanmean(data2(:,ii)); %calculate the column mean for each column
    dataRow = data2(:,ii); %pull out a row to currently work on
    dataRow(isnan(dataRow)) = columnAverage; %replace all NaN's in the column with the mean of that column
    cleanData(:,ii) = dataRow; %add revised row to new clean data matrix. 
end
handles.data = cleanData; %update the data handle to reflect the new clean data

studentAdmin = round(handles.data(:,1)); %caluclate the number of students addmitted from each school
totalAdmin = sum(handles.data(:,1)); %calculate the total number of students admitted form all schools
for ii = 1 : length(handles.data);
    handles.data(ii,1) = studentAdmin(ii)/totalAdmin ; %replace each number of students with the schools proportion of the total number of studnets admitted
end

set(handles.accStudSlide, 'Visible', 'On'); %turn on slider
set(handles.zeroPercent,'Visible','On'); %turn on left slider label
set(handles.twoPercent,'Visible','On'); %turn on center slider label
set(handles.fourPercent,'Visible','On'); %turn on right slider label
set(handles.sliderText,'Visible','On'); %turn on slider instructions
set(handles.recode, 'Visible','On'); %turn recode button on
set(handles.sliderVal, 'Visible','On'); %turn the slider output box on
set(handles.text16, 'Visible', 'On'); %turn reload instructions on

guidata(hObject, handles);


% --- Executes on button press in recode.
function recode_Callback(hObject, eventdata, handles)
% hObject    handle to recode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.PCA, 'Visible', 'On'); %turns PCA button on
set(handles.aboveTxt, 'Visible', 'On'); %turns recode output 'above' text boxes on
set(handles.belowTxt, 'Visible', 'On'); %turns recode output 'below' text boxes on
set(handles.numAbove, 'Visible', 'On'); %turns recode 'above' output boxes on
set(handles.numBelow, 'Visible', 'On'); %turns recode 'below' output boxes on
%histogram(handles.data(:,1))
above = length(find(handles.data(:,1) == 1)); %counts how many schools are a higher proportion of addmission than the slider 
below = length(find(handles.data(:,1) == 0)); %counts how many schools are a lower proportion of addmission than the slider
set(handles.numAbove, 'String', above); %shows how many schools are above slider position
set(handles.numBelow, 'String', below); %shows how many schools are below slider position

guidata(hObject, handles);


% --- Executes on slider movement.
function accStudSlide_Callback(hObject, eventdata, handles)
% hObject    handle to accStudSlide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
sliderInput = get(hObject, 'Value'); %gets the position of the slider
sliderInput2 = sliderInput*100; %multiplies the position location by 100 to match percentile labels
set(handles.sliderVal, 'String', sliderInput2); %shows the exact percent the slider is positioned on
set(handles.recode, 'Visible','On'); %turns recode button on
set(handles.sliderVal, 'Visible','On'); %turns the slider output box on
handles.data(handles.data(:,1) <= sliderInput) = 0; %updates all schools with admission proprotion lower than the slider to 1
handles.data(handles.data(:,1) > sliderInput) = 1; %updates all schools with admission proportion higher than the slider to 0
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function accStudSlide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to accStudSlide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'));
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

guidata(hObject, handles);


% --- Executes on button press in PCA.
function PCA_Callback(hObject, eventdata, handles)
% hObject    handle to PCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[loadings, dataNewCoor, eigVal] = pca(zscore(handles.data(:,2:16))); %makes a PCA from the z-scores of 15x predictors
handles.dataNewCoor = dataNewCoor;
eigenValMag = eig(corrcoef(handles.data)); %calculates the eigen value magnitudes 
bar(1:length(eigenValMag),sortrows(eigenValMag,-1)); %makes a bar graph of the eigen values (scree plot)
title ('Scree Plot');
line([min(xlim) max(xlim)], [1 1], 'linestyle', '--'); %kaiser criteron line
xlabel('Factors, in order of decreasing Eigenvalue');
ylabel('Eigenvalue');
box off;
set(handles.numFactors, 'Visible','On'); %turns the input box for # of factors on
set(handles.clustering, 'Visible', 'On'); %turns the clustering button on
set(handles.factTxt, 'Visible','On'); %turns the text about choosing factors on
set(handles.optBox, 'Visible','On'); %turns the check box on for optimal choice
set(handles.xGroup, 'Visible', 'On'); %turns the buttons on for the x group
set(handles.yGroup, 'Visible', 'On'); %turns the buttons on for the y group
set(handles.showFig, 'Visible', 'On'); %turns on the show figure button
set(handles.text13, 'Visible', 'On'); %turns on the instructions 
set(handles.clearFig, 'Visible', 'On'); %turn the clear figure button on
set(handles.classify, 'Visible', 'On'); %turn classify button on
set(handles.classButt, 'Visible', 'On'); %turn classify option buttons on
set(handles.text14, 'Visible', 'On'); %turn forest instructions on
set(handles.randForInput, 'Visible', 'On'); %turn # trees input on
set(handles.modelAcc, 'Visible', 'On'); %turn model accuracy output on
set(handles.valButt, 'Visible', 'On'); %turn validation option buttons on
guidata(hObject, handles);

% --- Executes on selection change in numFactors.
function numFactors_Callback(hObject, eventdata, handles)
% hObject    handle to numFactors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns numFactors contents as cell array
%        contents{get(hObject,'Value')} returns selected item from numFactors


% --- Executes during object creation, after setting all properties.
function numFactors_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numFactors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'));
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clustering.
function clustering_Callback(hObject, eventdata, handles)
% hObject    handle to clustering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[loadings, dataNewCoor, eigVal] = pca(zscore(handles.data(:,2:16))); %pull up PCA outputs

%get input for what loadings will be used in the clustering analysis. 
x1 = get(handles.xFact1, 'Value'); %if education is selected on x
x2 = get(handles.xFact2, 'Value'); % if community is selected on x
x3 = get(handles.xFact3, 'Value'); % if economy is selected on x
y1 = get(handles.yFact1, 'Value');%if education is selected on y
y2 = get(handles.yFact2, 'Value');% if community is selected on y
y3 = get(handles.yFact3, 'Value');% if economy is selected on y

if x1 == 1 %if education is selected on x
    xx = 1;  % pull from the 1st loading  
elseif x2 == 1 % if community is selected on x
    xx = 2;   %pull from the 2nd loading
elseif x3 == 1 % if economy is selected on x
    xx = 3; %pull fromt he 3rd loading
end

if y1 == 1 %if education is selected on y
    yy = 1; %pull from the 1st loading
elseif y2 == 1 % if community is selected on y
    yy = 2; %pull from the 2nd loading
elseif y3 == 1 % if economy is selected on y
    yy = 3; %pull from the 3rd loading
end

if get(handles.optBox, 'Value') == 0 %if optimal choice box is not selected
    
    numFacts = str2num(get(handles.numFactors,'String')); %use the number input into the box
    
elseif get(handles.optBox, 'Value') == 1 %if optimal choice box is selected
    
    %find the sillhouette scores
    factorX = [dataNewCoor(:,xx),dataNewCoor(:,yy)]; %make X the new data coordintes from selected factor loadings
    for ii = 2:16
        [cId cCoords SSd] = kmeans(factorX,ii); %find kmeans using factor chosen above
        s = silhouette(factorX,cId); %caluclate the sillhouette score
        sumSill(ii) = sum(s); %put summed sillhoette scores into matrix
        
    end
    %select the highest sillouhette score and what # of factors it corresponds to 
    maxSumSill = max(sumSill) %find the highest value of the matrix of summed sillhouette scores
    numFacts = find(sumSill==maxSumSill) %make the factor number the row of the highest sillhette score
end

factorX = [dataNewCoor(:,xx), dataNewCoor(:,yy)]; %make variable factorX the new coordinates of the factor loadings chosen

[cId, cCoords, SSd] = kmeans(factorX,numFacts); %run kmeans using the number of factors decided above
indexVector = [1:length(unique(cId))]; %create an index for loop

%plot all schools according to cluster color
for ii = indexVector
    plotIndex = find(cId == ii); %make index for plot
    plot(dataNewCoor(plotIndex,1),dataNewCoor(plotIndex,2),'.','markersize',10); %plot the data points for all schools
    hold on
    plot(cCoords(ii,1),cCoords(ii,2),'.','markersize',30,'color','k'); %Plot the cluster centers
end

guidata(hObject, handles);


% --- Executes on button press in optBox.
function optBox_Callback(hObject, eventdata, handles)
% hObject    handle to optBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of optBox
get(hObject,'Value')


% --- Executes on button press in classify.
function classify_Callback(hObject, eventdata, handles)
% hObject    handle to classify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 ll = length(handles.data); %find length of data
 predictors = handles.data(:,2:16); %isolates only the predictor variables
outcomes = handles.data(:,1); %isolates the outcome variable
if get(handles.val1, 'Value') == 1 %if the 80/20 option is selected
    randOrd = randperm(ll); %get a random order of numbers the length of the data matrix
    tt = round(ll*.8); %find what trial number is the 80% cut off
    ii = randOrd(1:tt); %80% of trials
    pp = randOrd(tt+1:ll); %20% of trials
    
    for rr = 1:tt %for trials 1 through the 80% cut off
    numSelect = ii(rr); %find the rr'th row of the randomized set of numbers (equal in length to rr)
    outTrain(rr,:) = outcomes(numSelect,:); %make a variable for the outputs to train on (80% of trials) by selecting them at random from the total trials
    predTrain(rr,:) = predictors(numSelect,:); %Make a matching variable of the predictors to use the same 80% of trials as the outputs
end
for rr = 1:length(pp) %for trials 1 through last 20%
    numSelect = pp(rr); %find the rr'th row of the randomized set of numbers (equal to rr in length)
    outTest(rr,:) = outcomes(numSelect,:); %create a variable of outputs to test on containing the other 20% of trials from the original trial set
    predTest(rr,:) = predictors(numSelect,:); %create a predictor variable to test on using the same 20% of trials 
end

elseif get(handles.val2, 'Value') == 1 %if the 'all but 1' option is selected
    randOrd = randperm(ll); %get a random order of numbers the length of the data matrix
    ii = randOrd(2:ll); %all trials except one (randomly decided)
    pp = randOrd(1); %one trial
   
    for rr = 2:ll-1 %two through lenth of L (minus 1- because starting on 2)
    numSelect = ii(rr); %find the rr'th row of the randomized set of numbers (equal to rr in length)
    outTrain(rr,:) = outcomes(numSelect,:); %create a variable for all but one of the trails to train predictions with
    predTrain(rr,:) = predictors(numSelect,:); %create a varaiable for the all but one of the trials to train outputs with
    end
    for rr = 1
    numSelect = ii(rr); %Select the first row of matrix ii
    outTest(rr,:) = outcomes(numSelect,:); %make a variable for the one test variable
    predTest(rr,:) = predictors(numSelect,:); %make a variable for the one train variable
end
end

 
    
if get(handles.svmButt, 'Value') == 1 %if the SVM model option is selected
    svmModel = fitcsvm(predTrain, outTrain); %creates the model using the amount of trials declared above
    [decision score] = predict(svmModel,predTest); %tests the model using the unseen predictors
    comp = [decision outTest];  %compares the model with the unseen outcomes vs decision
    modelAccuracy = (sum(comp(:,1) == comp(:,2))./length(comp)*100); %percentage of model accuracy
    num2str(set(handles.modelAcc, 'String', modelAccuracy)); %display the model accuracy
    
elseif get(handles.randFor, 'Value') == 1 %if the Random Forest option is selected
    numFor =  str2num(get(handles.randForInput, 'String')); %get the number of trees input by participant
    
    %run the treebagger model 
    treeModel = TreeBagger(numFor,predTrain,outTrain); %using the number of trials specified above
    [treeDecisions score] = predict(treeModel,predTest); %use the model to make a prediction using unseen predictors
    %The output of the treebagger are labels, we need to convert to numbers
    numericalOutcomes = str2num(char(string(treeDecisions))); %get the numerical output from the prediction
    empericalVsPrediction = [outTest numericalOutcomes]; %compare the unseen outcomes to the predicted values
    modelAccuracy = (sum(empericalVsPrediction(:,1) == empericalVsPrediction(:,2))./length(empericalVsPrediction)*100); %calculate percent model accuracy
    num2str(set(handles.modelAcc, 'String', modelAccuracy)); %display model accuracy
end
guidata(hObject, handles);

%%% Random Forest had the highest percent accuracy, especially when
%%% multiple trials are used to test the model. Since the output is
%%% categorical as 0 or 1, when only using 1x trail to predict the outcome
%%% there is a 50/50 chance of getting it correct. 

% --- Executes on button press in showFig.
function showFig_Callback(hObject, eventdata, handles)
% hObject    handle to showFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[loadings, dataNewCoor, eigVal] = pca(zscore(handles.data(:,2:16)));
x1 = get(handles.xFact1, 'Value'); %if education is selected on x
x2 = get(handles.xFact2, 'Value'); % if community is selected on x
x3 = get(handles.xFact3, 'Value'); % if economy is selected on x
y1 = get(handles.yFact1, 'Value');%if education is selected on y
y2 = get(handles.yFact2, 'Value');% if community is selected on y
y3 = get(handles.yFact3, 'Value');% if economy is selected on y

if x1 == 1 %if education is selected on x
    xx = 1;  % pull from the 1st loading  
elseif x2 == 1 % if community is selected on x
    xx = 2;   %pull from the 2nd loading
elseif x3 == 1 % if economy is selected on x
    xx = 3; %pull fromt he 3rd loading
end

if y1 == 1 %if education is selected on y
    yy = 1; %pull from the 1st loading
elseif y2 == 1 % if community is selected on y
    yy = 2; %pull from the 2nd loading
elseif y3 == 1 % if economy is selected on y
    yy = 3; %pull from the 3rd loading
end

plot(dataNewCoor(:,xx), dataNewCoor(:,yy), '.'); %make a scatter plot comparing the two loadings

guidata(hObject, handles);



function randForInput_Callback(hObject, eventdata, handles)
% hObject    handle to randForInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of randForInput as text
%        str2double(get(hObject,'String')) returns contents of randForInput as a double


% --- Executes during object creation, after setting all properties.
function randForInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to randForInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clearFig.
function clearFig_Callback(hObject, eventdata, handles)
% hObject    handle to clearFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.CA); % Make CA the current axes.
cla reset; %clear CA
