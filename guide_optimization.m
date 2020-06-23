function varargout = guide_optimization(varargin)
% GUIDE_OPTIMIZATION MATLAB code for guide_optimization.fig
%      GUIDE_OPTIMIZATION, by itself, creates a new GUIDE_OPTIMIZATION or raises the existing
%      singleton*.
%
%      H = GUIDE_OPTIMIZATION returns the handle to a new GUIDE_OPTIMIZATION or the handle to
%      the existing singleton*.
%
%      GUIDE_OPTIMIZATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIDE_OPTIMIZATION.M with the given input arguments.
%
%      GUIDE_OPTIMIZATION('Property','Value',...) creates a new GUIDE_OPTIMIZATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guide_optimization_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guide_optimization_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guide_optimization

% Last Modified by GUIDE v2.5 02-Jun-2020 01:05:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guide_optimization_OpeningFcn, ...
                   'gui_OutputFcn',  @guide_optimization_OutputFcn, ...
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


% --- Executes just before guide_optimization is made visible.
function guide_optimization_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guide_optimization (see VARARGIN)
clc;
VDD=0.2:0.2:0.6;
Wa=200*10^-9:200*10^-9:1500*10^-9;
[VDDI WAI] = meshgrid(VDD, Wa);
threshold = str2double(get(handles.edit7,'string'));

%threshold=0.49;
VT=(2.5.*threshold).*((85.*VDDI)./100);
%VT=(85.*VDDI)./100;
Wp=400*10^-9;
Vs=VDDI-VT;
Vgs=Vs; %%%%%%%%%%%%%%%may be changed%%%%%%%%%%%%%%%%
Wn=1000*10^-9;
L=60*10^-9;
%eta=0.1;
eta = str2double(get(handles.edit13,'string'));
Up=8.3*10^-3;
Un=18.9*10^-3;
Vcn=115*10^-3;
Ecn=Vcn/L;
Vdsat4=(36.*VDDI)./100;
Vds4=VDDI;
Vsatn=(Un*Ecn)/2;
VA=1;
r=Wn./WAI;
q=Wp./WAI;
ath=-VT;   %%%%%%%%%%%Vsb=0%%%%%%%%%%%%%%%
vt=0.026;
S=(vt^2)*exp(1.8);
C=1+(Vds4/VA);
Q=(C*L*Vsatn)/S;
x=exp(Vgs+ath);
y=exp(Vgs+ath+(VDDI*(1-eta)));
z=VDDI+VT-Vdsat4;
X=(Q.*z)./(S.*(((1./r).*x)+(q.*y)));
B=(log(X)./(log(exp(1)).*eta))-(log(r./q)./log(exp(1)))-(VDDI.*(1-eta));
Vr=exp(VDDI./VT).*(1./r);
a=-2*vt*B;
b=1-Vr+(2.*vt.*B.*(Vr+1));
c=2.*VT.*B.*r;
F=(b.^2)-(4.*a.*c);
T=(b-sqrt(F))./(a./4);
SNM=1.2.*vt.*(log(T)./log(exp(1)))+(Vs.*VT.*0.1.*(r.^2));
assignin('base','SNM',SNM);
clc;
N = str2double(get(handles.edit4,'string'));
%N=16;
%M=64;
M = str2double(get(handles.edit5,'string'));
qq= 1.602176565*10^-19;
KK= 1.380649*10^-23;
temp = str2double(get(handles.edit3,'String'));
assignin('base','KK',KK);
Vt=(KK.*(temp+273.15))./qq;
eox=3.45*10^-11;
%toxn=4*10^-9;
%toxp=4.2*10^-9;
toxn= str2double(get(handles.edit8,'string'));
toxp = str2double(get(handles.edit9,'string'));
Coxn=eox./toxn;
Coxp=eox./toxp;
Cjp=1.5*10^-3;
Cjn=2*10^-3; %1.27^10^-3%
Cjsn=1.27*10^-3;
Cjsp=1.06*10^-3;
nn=(Coxn+Cjn)./(Coxn);
np=(Coxp+Cjp)./(Coxp);
UP=8.6*10^-3;
Un=16.8*10^-3;
L=60*10^-9;
%Wn=200*10^-9:100*10^-9:800*10^-9;
Wn=800*10^-9;
Wp=200*10^-9;
Ya=0.091; %channel length modulation parameter%
%Yn=0.056; %channel length modulation parameter%
%Yp=2; %channel length modulation parameter%
Yn = str2double(get(handles.edit10,'string'));
Yp = str2double(get(handles.edit11,'string'));
Cjd=144*10^-19;
Vgsn=(85.*VDDI)./100;%
VTH=(83.*VDDI)./100;%
Vgsp=-(88.*VDDI)./100;%
Vgsa=436*10^-6;%
Vdsn=(5.*VDDI)./100;
Vdsp=(5.*VDDI)./100;
%eta=0.1; %drain induced barrier lowering%
eta = str2double(get(handles.edit13,'string'));
IDon=(Un.*Coxn.*(nn-1).*(Vt.^2).*(exp(-VTH./(nn.*Vt))).*(exp(-Vdsn./VTH)));
IDop=(UP.*Coxp.*(np-1).*(Vt.^2).*(exp(-VTH./(np.*Vt))).*(exp(-Vdsp./VTH)));
%Capacitance%
CBL=N.*Cjd;
Cgsn=(Wn.*L).*((Coxn.*Cjsn)./(Coxn+Cjsn));
Cgsa=(WAI.*L).*((Coxn.*Cjsn)./(Coxn+Cjsn));
Cgsp=(Wp.*L).*((Coxp.*Cjsp)./(Coxp+Cjsp));
CDI=M.*Cjd;  %M=number of columns 
CTX=CBL+Cgsa+Cgsp+Cgsn+CDI;
CTR=Cgsn+Cgsp+CBL+CDI;
%Subthreshold Currents%
IDn=IDon.*((Wn./L).*(exp(Vgsn/(nn.*Vt))));
IDp=IDop.*((Wp./L).*(exp(Vgsp/(np.*Vt))));
IDa=IDon.*((WAI./L).*(exp(Vgsn/(nn.*Vt))));
%Gm%
gmn=IDon.*(Wn./L).*(exp(Vgsn./(nn.*Vt))).*(1./(nn.*Vt));
gma=IDon.*(WAI./L).*(exp(Vgsn./(nn.*Vt))).*(1./(nn.*Vt));
gmp=IDop.*(Wp./L).*(exp(Vgsp./(np.*Vt))).*(1./(np.*Vt));
%Resistance%
ra=1./(Ya.*IDa);
rn=1./(Yn.*IDn);
rp=1./(Yp.*IDp);
RTX=(ra.*rn.*rp)./((ra.*(rn+rp))+(rn.*rp));
RTR=(rn.*rp)./(rn+rp);
%Read Access time for Array SRAM cell%
Vfinal=(15.*VDDI)./100;
Vintial=VDDI;
T=(log(Vintial./Vfinal)./log(exp(1)));
K=(1./(CTR.*RTR))-((gmn-gmp).*((1./CTR)+(1./CTX)))+((gma+(1./RTR))./CTX);
tawo=-1./K;
Delayr=tawo.*T;
clc;
Vs=VDDI+VT;
Vgs=Vs; 
Wp=400*10^-9;
L=60*10^-9;
Up=8.6*10^-3;
Un=18.9*10^-3;
Vcn=115*10^-3;
Vcp=115*10^-3;
Ecn=Vcn/L;
Vdsat3=(25.*VDDI)./100;
Vdsat4=(20.*VDDI)./100;
Vds3=VDDI;
Vds4=VDDI;
Vsatn=(Un*Ecn)/2;
VA3=0.88;
VA4=0.88;
B5=(Up*Vcp)/L;
B3=Vsatn*(1+(Vds3/VA3));
B4=Vsatn*(1+(Vds4/VA4));
B2=(Un*Vcn)/L;
r=Wn./WAI;
q=Wp./WAI;
a=B4./(r.*B2);
d=B3./(q.*B5);
%%%%%%%%%%%%%%%%Equation%%%%%%%%%%%%%%%%%%%%%
Vr=VT-(a.*Vs)+(a.*Vdsat4)+(a.*Vcn);
VQ=4.*a.*Vcn.*(a-0.5).*(Vs-Vdsat4);
S=(Vs+Vr).^2+VQ;
A=(2.*(Vs+Vr))./(sqrt(S));
K=(-1./((2.*a)+1)).*(-1-A);  %%%%%%%%%%%%%%%may be changed%%%%%%%%%%%%%%%%
Y=-(Vs+Vr)+sqrt(S);
V0=(2.*Y./(2.*(a)))+(K.*Vs);  %%%%%%%%%%%%%%%may be changed%%%%%%%%%%%%%%%%
WNM=Vdsat3.*(VDDI+VT-(K.*Vgs)+V0+(0.1./q.*Vs)+(((0.1./q)).*Vdsat3)+(0.2.*sqrt(1./q.^2)));
%%%%%%%%%%%%%%%%%%%%%Plotting%%%%%%%%%%%%%%%%%%%%

clc;
%VT=0.36;
%N=16;
N = str2double(get(handles.edit4,'string'));
%M=64;
M = str2double(get(handles.edit5,'string'));
KK= 1.380649*10^-23;
temp = str2double(get(handles.edit3,'String'));
qq= 1.602176565*10^-19;
Vt=(KK.*(temp+273.15))./qq;
%Vt=0.026;
%VDDI=0.5;
eox=3.45*10^-11;
%toxn=4*10^-9;
%toxp=4.2*10^-9;
toxp = str2double(get(handles.edit9,'string'));
toxn= str2double(get(handles.edit8,'string'));
Coxn=eox./toxn;
Coxp=eox./toxp;
Cjp=1.5*10^-3;
Cjn=2*10^-3; %1.27^10^-3%
Cjsn=1.27*10^-3;
Cjsp=1.06*10^-3;
nn=(Coxn+Cjn)./(Coxn);
np=(Coxp+Cjp)./(Coxp);
Up=8.6*10^-3;
Un=16.8*10^-3;
L=60*10^-9;
%Wn=200*10^-9:100*10^-9:800*10^-9;
Wp=400*10^-9;
Ya=0.091; %channel length modulation parameter%
%Yn=0.056; %channel length modulation parameter%
%Yp=2; %channel length modulation parameter%
Yp = str2double(get(handles.edit11,'string'));
Yn = str2double(get(handles.edit10,'string'));
Cjdn=140*10^-19;
Cjdm=120*10^-19;
VTHn=(85.*VDDI)./100;%
%VTHp=-(85.*VDDI)./100;%
thresholdp = str2double(get(handles.edit12,'string'));

VTHp=-(2.5.*thresholdp).*((85.*VDDI)./100);
Vgsn=(81.*VDDI)./100;%
Vgsp=-(81.*VDDI)./100;%
Vgsa=(81.*VDDI)./100;
%Vgsa=436*10^-6;%
Vdsn=(95.*VDDI)./100;%%no change in delay Vs SRAM ARRAY%
Vdsp=(-95.*VDDI)./100;%%no change in delay Vs SRAM ARRAY%
etaa=0.15; %drain induced barrier lowering%
%Capacitance%
CBL=N.*Cjdn;
Cgsn=(Wn.*L).*((Coxn.*Cjsn)./(Coxn+Cjsn));
Cgsa=(WAI.*L).*((Coxn.*Cjsn)./(Coxn+Cjsn));
Cgsp=(Wp.*L).*((Coxp.*Cjsp)./(Coxp+Cjsp));
CDI=M.*Cjdm;  %M=number of columns 
c=Cgsp+Cgsn;
Ca=CDI+CBL+Cgsa;
%Subthreshold Currents%
IDon=Un.*Coxn.*(nn-1).*(Wn./L).*(Vt.^2).*(exp((Vgsn-VTHn)./(nn.*Vt)));
IDop=Up.*Coxp.*(np-1).*(Wp./L).*(Vt.^2).*(exp((Vgsp-VTHp)./(np.*Vt)));
IDoa=Un.*Coxn.*(nn-1).*(WAI./L).*(Vt.^2).*(exp((Vgsa-VTHn)./(nn.*Vt)));
IDn=IDon.*(1-(exp(-Vdsn./Vt))).*(exp((etaa.*Vdsn)./(nn.*Vt)));%
IDp=IDop.*(1-(exp(-Vdsp./Vt))).*(exp((etaa.*Vdsp)./(np.*Vt)));%
IDa=IDoa.*(1-(exp(-Vdsn./Vt))).*(exp((etaa.*Vdsn)./(nn.*Vt)));%
%%%%%%Transconductance-Gm%%%%%%%%
%gmn=IDon.*(Wn./L).*(exp(Vgsn./(nn.*Vt))).*(1./(nn.*Vt));
%gma=IDoa.*(Wa./L).*(exp(Vgsn./(nn.*Vt))).*(1./(nn.*Vt));
%gmp=IDop.*(popupmenu1./L).*(exp(Vgsp./(np.*Vt))).*(1./(np.*Vt));
gmn=IDn./(nn.*Vt); %%Micheal Perrot slides%%%
gmp=IDp./(np.*Vt);
gma=IDa./(nn.*Vt);
X=gmn-gmp;
%%%%%%%%%%%%%Resistance%%%%%%%%%%%%
%ra6=1./(Ya.*IDa);%1
%ra6=VDD./IDa;%2
A=exp((etaa.*VDDI)./(2.*nn.*Vt));%
B=exp((etaa.*VDDI)./(nn.*Vt));%
AA=((etaa.*VDDI)./(2))+(nn.*Vt);%
BB=(etaa.*VDDI)+(nn.*Vt);%
Const=(2.*nn.*Vt)./(IDa.*((etaa)^2).*VDDI);
ra6=Const.*((AA./A)-(BB./B));
rn=1./(Yn.*IDn);
rp=1./(Yp.*IDp);
r=(rn.*rp)./(rn+rp);
%Read Access time for Array SRAM cell%
Vfinal=(25.*VDDI)./100;
%Vfinal=(91.*VDDI)./100;
Vintial=VDDI;
T=(log(Vintial./Vfinal)./log(exp(1)));
ln2=(log(2)./log(exp(1)));
K=(X-(1./r))./(c);
tawo1=1./K;
tawo2=(ra6.*Ca);
tawo=tawo1+tawo2;
Delayw=tawo.*T;
Delay=0.5.*(Delayr+Delayw);
NM=0.5.*(SNM+WNM);
FOM=Delay./NM;
axes(handles.axes2)
handles.output = hObject;
handles.vdd= VDDI;
handles.Wa=WAI;
handles.FOM=FOM;
   
% Choose default command line output for guide_optimization
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes guide_optimization wait for user response (see UIRESUME)
% uiwait(handles.figure1);
%uiwait(fig);
% --- Outputs from this function are returned to the command line.
function varargout = guide_optimization_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%UIRESUME(handles.guide_optimization);
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
VDD=0.2:0.2:0.6;
Wa=200*10^-9:200*10^-9:1500*10^-9;
[VDDI WAI] = meshgrid(VDD, Wa);
%threshold=0.49;
threshold = str2double(get(handles.edit7,'string'));
assignin('base','threshold',threshold);
VT=(2.5.*threshold).*((85.*VDDI)./100);
%VT=(85.*VDDI)./100;
N = str2double(get(handles.edit4,'string'));
assignin('base','N',N);
%N=16;
M = str2double(get(handles.edit5,'string'));
assignin('base','M',M);
%M=64;
eox=3.45*10^-11;
%toxn=4*10^-9;
%toxp=4.2*10^-9;
toxp = str2double(get(handles.edit9,'string'));
toxn= str2double(get(handles.edit8,'string'));
Coxn=eox./toxn;
assignin('base','Coxn',Coxn);
Coxp=eox./toxp;
assignin('base','Coxp',Coxp);
Cjp=1.5*10^-3;
Cjn=2*10^-3; %1.27^10^-3%
Cjsn=1.27*10^-3;
Cjsp=1.06*10^-3;
nn=(Coxn+Cjn)./(Coxn);
np=(Coxp+Cjp)./(Coxp);
UP=8.6*10^-3;
Un=16.8*10^-3;
L=60*10^-9;
%Wn=200*10^-9:100*10^-9:800*10^-9;
Wn=800*10^-9;
Wp=200*10^-9;
Ya=0.091; %channel length modulation parameter%
%Yn=0.056; %channel length modulation parameter%
%Yp=2; %channel length modulation parameter%
Yp = str2double(get(handles.edit11,'string'));
Yn = str2double(get(handles.edit10,'string'));
assignin('base','Yp',Yp);
assignin('base','Yn',Yn);
Cjd=144*10^-19;
Vgsn=(85.*VDDI)./100;%
VTH=(83.*VDDI)./100;%
Vgsp=-(88.*VDDI)./100;%
Vgsa=436*10^-6;%
Vdsn=(5.*VDDI)./100;
Vdsp=(5.*VDDI)./100;
%eta=0.1; %drain induced barrier lowering%
eta = str2double(get(handles.edit13,'string'));
temp = str2double(get(handles.edit3,'string'));
assignin('base','temp',temp);
qq= 1.602176565*10^-19;
KK= 1.380649*10^-23;
Vt= (KK.*(temp+273.15))./qq;
xxx= Vt;
IDon=(Un.*Coxn.*(nn-1).*(Vt.^2).*(exp(-VTH./(nn.*xxx))).*(exp(-Vdsn./VTH)));
IDop=(UP.*Coxp.*(np-1).*(Vt.^2).*(exp(-VTH./(np.*xxx))).*(exp(-Vdsp./VTH)));

assignin('base','Vt',Vt);
assignin('base','xxx',xxx);
assignin('base','IDop',IDop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Capacitance%
CBL=N.*Cjd;
Cgsn=(Wn.*L).*((Coxn.*Cjsn)./(Coxn+Cjsn));
Cgsa=(WAI.*L).*((Coxn.*Cjsn)./(Coxn+Cjsn));
Cgsp=(Wp.*L).*((Coxp.*Cjsp)./(Coxp+Cjsp));
CDI=M.*Cjd;  %M=number of columns 
CTX=CBL+Cgsa+Cgsp+Cgsn+CDI;
CTR=Cgsn+Cgsp+CBL+CDI;
%Subthreshold Currents%
IDn=IDon.*((Wn./L).*(exp(Vgsn/(nn.*Vt))));
IDp=IDop.*((Wp./L).*(exp(Vgsp/(np.*Vt))));
IDa=IDon.*((WAI./L).*(exp(Vgsn/(nn.*Vt))));
assignin('base','IDa',IDa);
%Gm%
gmn=IDon.*(Wn./L).*(exp(Vgsn./(nn.*Vt))).*(1./(nn.*Vt));
gma=IDon.*(WAI./L).*(exp(Vgsn./(nn.*Vt))).*(1./(nn.*Vt));
gmp=IDop.*(Wp./L).*(exp(Vgsp./(np.*Vt))).*(1./(np.*Vt));
assignin('base','gmp',gmp);
%Resistance%
ra=1./(Ya.*IDa);
rn=1./(Yn.*IDn);
rp=1./(Yp.*IDp);
RTX=(ra.*rn.*rp)./((ra.*(rn+rp))+(rn.*rp));
assignin('base','RTX',RTX);
RTR=(rn.*rp)./(rn+rp);
%Read Access time for Array SRAM cell%
Vfinal=(15.*VDDI)./100;
Vintial=VDDI;
T=(log(Vintial./Vfinal)./log(exp(1)));
assignin('base','T',T);
K=(1./(CTR.*RTR))-((gmn-gmp).*((1./CTR)+(1./CTX)))+((gma+(1./RTR))./CTX);
assignin('base','K',K);
tawo=-1./K;
Delayr=tawo.*T;
assignin('base','Delayr',Delayr);
Cjdn=140*10^-19;
Cjdm=120*10^-19;
VTHn=(85.*VDDI)./100;%
%VTHp=-(85.*VDDI)./100;%
thresholdp = str2double(get(handles.edit12,'string'));
assignin('base','thresholdp',thresholdp);
VTHp=-(2.5.*thresholdp).*((85.*VDDI)./100);
Vgsn=(81.*VDDI)./100;%
Vgsp=-(81.*VDDI)./100;%
Vgsa=(81.*VDDI)./100;
Vdsn=(95.*VDDI)./100;%%no change in delay Vs SRAM ARRAY%
Vdsp=(-95.*VDDI)./100;%%no change in delay Vs SRAM ARRAY%
etaa=0.15; %drain induced barrier lowering%
%Capacitance%
CBL=N.*Cjdn;
Cgsn=(Wn.*L).*((Coxn.*Cjsn)./(Coxn+Cjsn));
Cgsa=(WAI.*L).*((Coxn.*Cjsn)./(Coxn+Cjsn));
Cgsp=(Wp.*L).*((Coxp.*Cjsp)./(Coxp+Cjsp));
CDI=M.*Cjdm;  %M=number of columns 
c=Cgsp+Cgsn;
Ca=CDI+CBL+Cgsa;
%Subthreshold Currents%
IDon=Un.*Coxn.*(nn-1).*(Wn./L).*(Vt.^2).*(exp((Vgsn-VTHn)./(nn.*Vt)));
IDop=UP.*Coxp.*(np-1).*(Wp./L).*(Vt.^2).*(exp((Vgsp-VTHp)./(np.*Vt)));
IDoa=Un.*Coxn.*(nn-1).*(WAI./L).*(Vt.^2).*(exp((Vgsa-VTHn)./(nn.*Vt)));
assignin('base','IDoa',IDoa);
IDn=IDon.*(1-(exp(-Vdsn./Vt))).*(exp((etaa.*Vdsn)./(nn.*Vt)));%
IDp=IDop.*(1-(exp(-Vdsp./Vt))).*(exp((etaa.*Vdsp)./(np.*Vt)));%
IDa=IDoa.*(1-(exp(-Vdsn./Vt))).*(exp((etaa.*Vdsn)./(nn.*Vt)));%
%%%%%%Transconductance-Gm%%%%%%%%
gmn=IDn./(nn.*Vt); %%Micheal Perrot slides%%%
gmp=IDp./(np.*Vt);
gma=IDa./(nn.*Vt);
X=gmn-gmp;
%%%%%%%%%%%%%Resistance%%%%%%%%%%%%
A=exp((etaa.*VDDI)./(2.*nn.*Vt));%
B=exp((etaa.*VDDI)./(nn.*Vt));%
AA=((etaa.*VDDI)./(2))+(nn.*Vt);%
BB=(etaa.*VDDI)+(nn.*Vt);%
Const=(2.*nn.*Vt)./(IDa.*((etaa)^2).*VDDI);
ra6=Const.*((AA./A)-(BB./B));
rn=1./(Yn.*IDn);
rp=1./(Yp.*IDp);
r=(rn.*rp)./(rn+rp);
%Read Access time for Array SRAM cell%
Vfinal=(25.*VDDI)./100;
Vintial=VDDI;
T=(log(Vintial./Vfinal)./log(exp(1)));
ln2=(log(2)./log(exp(1)));
K=(X-(1./r))./(c);
tawo1=1./K;
tawo2=(ra6.*Ca);
tawo=tawo1+tawo2;
Delayw=tawo.*T;
assignin('base','Delayw',Delayw);
Delay=0.5.*(Delayr+Delayw);
%%%%%%%%%%%%
Vs=VDDI+VT;
Vgs=Vs; 
Wp=400*10^-9;
L=60*10^-9;
Up=8.6*10^-3;
Un=18.9*10^-3;
Vcn=115*10^-3;
Vcp=115*10^-3;
Ecn=Vcn/L;
Vdsat3=(25.*VDDI)./100;
Vdsat4=(20.*VDDI)./100;
Vds3=VDDI;
Vds4=VDDI;
Vsatn=(Un*Ecn)/2;
VA3=0.88;
VA4=0.88;
B5=(Up*Vcp)/L;
B3=Vsatn*(1+(Vds3/VA3));
B4=Vsatn*(1+(Vds4/VA4));
B2=(Un*Vcn)/L;
r=Wn./WAI;
q=Wp./WAI;
a=B4./(r.*B2);
d=B3./(q.*B5);
%%%%%%%%%%%%%%%%Equation%%%%%%%%%%%%%%%%%%%%%
Vr=VT-(a.*Vs)+(a.*Vdsat4)+(a.*Vcn);
VQ=4.*a.*Vcn.*(a-0.5).*(Vs-Vdsat4);
S=(Vs+Vr).^2+VQ;
A=(2.*(Vs+Vr))./(sqrt(S));
K=(-1./((2.*a)+1)).*(-1-A);  %%%%%%%%%%%%%%%may be changed%%%%%%%%%%%%%%%%
Y=-(Vs+Vr)+sqrt(S);
V0=(2.*Y./(2.*(a)))+(K.*Vs);  %%%%%%%%%%%%%%%may be changed%%%%%%%%%%%%%%%%
WNM=Vdsat3.*(VDDI+VT-(K.*Vgs)+V0+(0.1./q.*Vs)+(((0.1./q)).*Vdsat3)+(0.2.*sqrt(1./q.^2)));
%%%%%%%%%%%%

Vs=VDDI-VT;
Vgs=Vs; %%%%%%%%%%%%%%%may be changed%%%%%%%%%%%%%%%%
Wn=1000*10^-9;
L=60*10^-9;
%eta=0.1;
eta = str2double(get(handles.edit13,'string'));
Up=8.3*10^-3;
Un=18.9*10^-3;
Vcn=115*10^-3;
Ecn=Vcn/L;
Vdsat4=(36.*VDDI)./100;
Vds4=VDDI;
Vsatn=(Un*Ecn)/2;
VA=1;
r=Wn./WAI;
q=Wp./WAI;
ath=-VT;   %%%%%%%%%%%Vsb=0%%%%%%%%%%%%%%%
vt=0.026;
S=(vt^2)*exp(1.8);
C=1+(Vds4/VA);
Q=(C*L*Vsatn)/S;
x=exp(Vgs+ath);
y=exp(Vgs+ath+(VDDI*(1-eta)));
z=VDDI+VT-Vdsat4;
X=(Q.*z)./(S.*(((1./r).*x)+(q.*y)));
B=(log(X)./(log(exp(1)).*eta))-(log(r./q)./log(exp(1)))-(VDDI.*(1-eta));
Vr=exp(VDDI./VT).*(1./r);
a=-2*vt*B;
b=1-Vr+(2.*vt.*B.*(Vr+1));
c=2.*VT.*B.*r;
F=(b.^2)-(4.*a.*c);
T=(b-sqrt(F))./(a./4);
SNM=1.2.*vt.*(log(T)./log(exp(1)))+(Vs.*VT.*0.1.*(r.^2));
%%%%%%%%%%%%%%%
NM=0.5.*(SNM+WNM);
FOM=(Delay./NM).*10^9;
assignin('base','FOM',FOM);
axes(handles.axes2)
handles.output = hObject;
handles.vdd= VDDI;
handles.Wa=WAI*10^6;
handles.FOM=FOM;
fig=surf(handles.vdd, handles.Wa, handles.FOM);
%title('Optimization FOM','Fontsize',50,'color','blue')
set (gca , 'fontsize', 9) ;
set (gca , 'fontweight', 'bold' ) ;
set (gca , 'linewidth', 2) ;
xlabel('\bf\fontname{helvetica}\fontsize{10}VDD (V)');
ylabel('\bf\fontname{helvetica}\fontsize{10}W_ access (\mum)');
zlabel('\bf\fontname{helvetica}\fontsize{10}FOM (ns/m)');
grid on;
box on;
pbaspect([2 1 2]);
guidata(hObject,handles);


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
%temp = get(handles.edit3, 'Value');
temp = str2double(get(handles.edit3,'string'));
%assignin('base','temp',temp);

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
N = str2double(get(handles.edit4,'string'));

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
M = str2double(get(handles.edit5,'string'));

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
threshold = str2double(get(handles.edit7,'string'));

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
toxn = str2double(get(handles.edit8,'string'));

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double
toxp = str2double(get(handles.edit9,'string'));

% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function text13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double
Yn = str2double(get(handles.edit10,'string'));

% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double
Yp = str2double(get(handles.edit11,'string'));

% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double
thresholdp = str2double(get(handles.edit12,'string'));

% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double
eta = str2double(get(handles.edit13,'string'));

% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
