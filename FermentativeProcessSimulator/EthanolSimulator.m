function varargout = EthanolSimulator(varargin)
% ETHANOLSIMULATOR MATLAB code for EthanolSimulator.fig
%      ETHANOLSIMULATOR, by itself, creates a new ETHANOLSIMULATOR or raises the existing
%      singleton*.
%
%      H = ETHANOLSIMULATOR returns the handle to a new ETHANOLSIMULATOR or the handle to
%      the existing singleton*.
%
%      ETHANOLSIMULATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ETHANOLSIMULATOR.M with the given input arguments.
%
%      ETHANOLSIMULATOR('Property','Value',...) creates a new ETHANOLSIMULATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EthanolSimulator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EthanolSimulator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EthanolSimulator

% Last Modified by GUIDE v2.5 14-May-2017 12:55:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EthanolSimulator_OpeningFcn, ...
                   'gui_OutputFcn',  @EthanolSimulator_OutputFcn, ...
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


% --- Executes just before EthanolSimulator is made visible.
function EthanolSimulator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EthanolSimulator (see VARARGIN)

% Choose default command line output for EthanolSimulator
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EthanolSimulator wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EthanolSimulator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function input1_editText_Callback(hObject, eventdata, handles)
% hObject    handle to input1_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
input = (get(hObject,'String'));
% Hints: get(hObject,'String') returns contents of input1_editText as text
%        str2double(get(hObject,'String')) returns contents of input1_editText as a double


% --- Executes during object creation, after setting all properties.
function input1_editText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input1_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Calcular_Button.
function Calcular_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Calcular_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
        global esc;
        %guidata(hObject,handles);
		%esc = str2num(get(handles.Modelo_editText, 'String'));
		
        global T;
        global X;
        global S;
        global P;
        global X1;
        global S1;
        global P1;
        X0 = X(1);
		S0 = S(1);
		P0 = 0;
				
		% Fatores de Conversão
		YSX = (S(1) - S(end))/(X(end) - X(1));
		YPX = (P(end) - P(1))/(X(end) - X(1));
				
		% Interpolação de Splines
        polyX = @(t) interp1(T, X, t, 'spline');
        polyS = @(t) interp1(T, S, t, 'spline');
        polyP = @(t) interp1(T, P, t, 'spline');
                
		% Valores máximos de Concentração
        [~, Xmax] = fminbnd(@(t) -polyX(t), T(1), T(end));
        Xmax = -Xmax;
               
        [~, Smax] = fminbnd(@(t) -polyS(t), T(1), T(end));
        Smax = -Smax;
                
        [~, Pmax] = fminbnd(@(t) -polyP(t), T(1), T(end));
        Pmax = -Pmax;
                
		%% Calculo dos Parâmetros dos Modelos
				
		% Discretização do tempo descartando os tempos cuja concentração de Substrato é nula
		h     = 0.001;
		t     = T(1):h:T(min(find(S,1,'last') + 1, length(T)));
				
		% Valores de concentração de Biomassa no Tempo discretizado
		X_mu  = polyX(t);
				
		% Aproximação da Derivada por diferenças finitas
		dX_mu = diff(X_mu)/h;
				
		% Calculo da Velocidade Específica de Transformação da Biomassa (mu_X)
		mu_X  = dX_mu./X_mu(1:end-1);
		t     = t(1:end-1);
		X_mu  = X_mu(1:end-1);
				
		% Valores da concentração de Substrato no Tempo
		S_mu  = polyS(t);
				
		% Determinação do valor maximo de mu_X
		[MAX_mu_X, ~] = max(mu_X);
		aux = MAX_mu_X/2; % Metade do Valor Máximo de mu_X
				
		% Determinar a posição de aux (Valor Máximo dividido por 2)
		for i=1:length(mu_X)
			if abs(mu_X(i) - aux) <= 1e-4
				pos = i;
			end
		end
				
		% Constante de Saturação
		KS = S_mu(pos);
				
		% Constante de Inibição pelo Produto e parâmetro n
		Ki = 1;
		n  = 1;
		
		polySt2   = polyS(t);
		polyPt2   = polyP(t);
		
		t = cputime;
        % Otimização dos Parâmetros				
		coef = lsqnonlin(@(miks_kin)calc_coef(miks_kin, mu_X, X, polySt2, polyPt2, Pmax, esc),[MAX_mu_X, KS, Ki, n]);
				
        x = lsqnonlin(@(miks)Otimizador(miks, T, X, S, P, t, Pmax, YSX, YPX, esc), coef);
						
		mu_MAX = x(1);
		KS     = x(2);
		Ki     = x(3);
		n      = x(4);
						
		% Resolução do Sistema de EDO's com os melhores parâmetros possíveis.
		[t1, X1, S1, P1] = RungeKutta4th('Biomassa', 'Substrato', 'Produto', T(1), T(end), 0.001, X0, S0, P0, mu_MAX, KS, Ki, n, Pmax, YSX, YPX, esc);
		
		time = cputime-t;
		
        set(handles.Time_staticText,'String',time);
        
		X1 = X1';
		S1 = S1';
		P1 = P1';
		auxX = X1;
		auxS = S1;
		auxP = P1;
		X1 = X1(1:2000:T(end)*1000+1);
		S1 = S1(1:2000:T(end)*1000+1);
		P1 = P1(1:2000:T(end)*1000+1);
        
        % Cálculo do Desvio Padrão Residual
		dpX = RSD(X, X1, length(X));
		dpS = RSD(S, S1, length(S));
		dpP = RSD(P, P1, length(P));
        
        % Teste de Confiança utilizando nível de Randomicidade
		Z1 = Randomic(X, X1);
		Z2 = Randomic(S, S1);
		Z3 = Randomic(P, P1);
        
        if(esc == 1)
            set(handles.Mu_staticText,'String',mu_MAX);
            set(handles.Ks_staticText,'String',KS);
        end
        if(esc == 2)
            set(handles.Mu_staticText,'String',mu_MAX);
            set(handles.Ks_staticText,'String',KS);
            set(handles.Ki_staticText,'String',Ki);
        end
        if(esc == 3)
            set(handles.Mu_staticText,'String',mu_MAX);
            set(handles.Ks_staticText,'String',KS);
            set(handles.n_staticText,'String',n);
        end
        if(esc == 4)
            set(handles.Mu_staticText,'String',mu_MAX);
            set(handles.Ks_staticText,'String',KS);
            set(handles.Ki_staticText,'String',Ki);
        end
        if(esc == 5)
            set(handles.Mu_staticText,'String',mu_MAX);
            set(handles.Ks_staticText,'String',KS);
        end
        
        set(handles.randX_staticText,'String',Z1);
        set(handles.randS_staticText,'String',Z2);
        set(handles.randP_staticText,'String',Z3);
        set(handles.rsdX_staticText,'String',dpX);
        set(handles.rsdS_staticText,'String',dpS);
        set(handles.rsdP_staticText,'String',dpP);
        
        guidata(hObject,handles);

function Modelo_editText_Callback(hObject, eventdata, handles)
% hObject    handle to Modelo_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
input = (get(hObject,'String'));
% Hints: get(hObject,'String') returns contents of Modelo_editText as text
%        str2double(get(hObject,'String')) returns contents of Modelo_editText as a double


% --- Executes during object creation, after setting all properties.
function Modelo_editText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Modelo_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Monod_pushbutton.
function Monod_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Monod_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global esc;
esc = 1;
guidata(hObject,handles);


% --- Executes on button press in Andrews_pushbutton.
function Andrews_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Andrews_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global esc;
esc = 2;
guidata(hObject, handles);

% --- Executes on button press in Levesnpiel_pushbutton.
function Levesnpiel_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Levesnpiel_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global esc;
esc = 3;
guidata(hObject, handles);

% --- Executes on button press in GhoseThyagi_pushbutton.
function GhoseThyagi_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to GhoseThyagi_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global esc;
esc = 4;
guidata(hObject, handles);

% --- Executes on button press in Jin_pushbutton.
function Jin_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Jin_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global esc;
esc = 5;
guidata(hObject, handles);

% --- Executes on button press in ok_pushbutton.
function ok_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ok_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
        global T;
        global X;
        global S;
        global P;
        a = get(handles.input1_editText, 'String');
		arquivo = fopen(a);
		A = fscanf(arquivo, '%f');
		fclose(arquivo);
		T = A(1:7);
		X = A(8:14);
		S = A(15:21);
		P = A(22:28);
        guidata(hObject,handles);


% --- Executes on button press in plot_pushbutton.
function plot_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to plot_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
		% Plotagem dos Dados
		global X;
		global S;
		global P;
		global T;
		global X1;
		global P1;
		global S1;
		fig = figure(); 
		set(fig,'Position',[240 240 900 400]);
		
		subplot(1,3,1); 
        hold on;  
		plot(T, X1, '--b', 'LineWidth', 1);
		plot(T, X, 'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b');
		title('Biomassa');
		ylabel('Concentração de Biomassa (g/L)');
		xlabel('Tempo (h)');
		legend('X (Numérico)', 'X (Empírico)', 'Location', 'southeast')

		subplot(1,3,2); 
        hold on;  
		plot(T, S1, '--k', 'LineWidth', 1);
		plot(T, S, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
		title('Substrato');
		ylabel('Concentração de Substrato (g/L)');
		xlabel('Tempo (h)');
		legend('S (Numérico)', 'S (Empírico)', 'Location', 'northeast')
		
		subplot(1,3,3); 
        hold on;  
		plot(T, P1, '--r', 'LineWidth', 1);
		plot(T, P, '^', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
		title('Produto');
		ylabel('Concentração de Produto (g/L)');
		xlabel('Tempo (h)');
		legend('P (Numérico)', 'P (Empírico)', 'Location', 'southeast')
		
		drawnow
		guidata(hObject,handles);
