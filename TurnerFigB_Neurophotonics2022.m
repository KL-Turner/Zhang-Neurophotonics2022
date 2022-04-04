%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate figure panel 3 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

colorRfcAwake = [(0/256),(64/256),(64/256)];
colorRfcNREM = [(0/256),(174/256),(239/256)];
colorRfcREM = [(190/256),(30/256),(45/256)];
colorOrange = [(233/256),(105/256),(44/256)];
colorPink = [(256/256),(28/256),(207/256)];
colorSapphire = [(15/256),(82/256),(187/256)];
colorDarkRed = [(164/256),(0/256),(0/256)];
exampleMergedFileID = 'T118_RH_191211_14_44_20_015_A1_MergedData.mat';
load(exampleMergedFileID,'-mat')
exampleSpecFileID = 'T118_RH_191211_14_44_20_015_A1_SpecData.mat';
load(exampleSpecFileID,'-mat')
exampleBaselinesFileID = 'T118_RestingBaselines.mat';
load(exampleBaselinesFileID,'-mat')
vesselID = 'A1';
strDay = 'Dec11';
dsFs = MergedData.notes.dsFs;
p2Fs = MergedData.notes.p2Fs;
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1,k1] = butter(4,10/(MergedData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
[z2,p2,k2] = butter(4,0.5/(MergedData.notes.p2Fs/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);
% whisker angle
filtWhiskerAngle = filtfilt(sos1,g1,MergedData.data.whiskerAngle);
% pressure sensor
filtForceSensor = filtfilt(sos1,g1,abs(MergedData.data.forceSensorL));
% EMG
EMG = MergedData.data.EMG.data;
normEMG = EMG - RestingBaselines.manualSelection.EMG.data.(strDay);
filtEMG = filtfilt(sos1,g1,normEMG);
% vessel diameter
vesselDiameter = MergedData.data.vesselDiameter.data;
normVesselDiameter = (vesselDiameter - RestingBaselines.manualSelection.vesselDiameter.data.(vesselID).(strDay))./(RestingBaselines.manualSelection.vesselDiameter.data.(vesselID).(strDay));
filtVesselDiameter = filtfilt(sos2,g2,normVesselDiameter)*100;
% cortical and hippocampal spectrograms
cortNormS = SpecData.corticalNeural.fiveSec.normS.*100;
hipNormS = SpecData.hippocampalNeural.fiveSec.normS.*100;
T = SpecData.corticalNeural.fiveSec.T;
F = SpecData.corticalNeural.fiveSec.F;
%% 2PLSM sleep example
figure;
% EMG and force sensor
ax1 = subplot(6,1,1);
p1 = plot((1:length(filtEMG))/dsFs,filtEMG,'color','k','LineWidth',0.5);
ylabel({'EMG','power (a.u.)'})
ylim([-2,2.5])
yyaxis right
p2 = plot((1:length(filtForceSensor))/dsFs,filtForceSensor,'color',colorPink,'LineWidth',0.5);
ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p1,p2],'EMG','Pressure')
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([15,75,135,195,255,315,375,435,495,555,615])
xlim([15,615])
ylim([-0.1,2.5])
ax1.TickLength = [0.01,0.01];
ax1.YAxis(1).Color = 'k';
ax1.YAxis(2).Color = colorPink;
% whisker angle
ax2 = subplot(6,1,2);
plot((1:length(filtWhiskerAngle))/dsFs,-filtWhiskerAngle,'color','k','LineWidth',0.5)
ylabel({'Whisker','angle (deg)'})
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([15,75,135,195,255,315,375,435,495,555,615])
ax2.TickLength = [0.01,0.01];
xlim([15,615])
ylim([-10,40])
% vessel diameter
ax34 = subplot(6,1,[3,4]);
p3 = plot((1:length(filtVesselDiameter))/p2Fs,filtVesselDiameter,'color',colorDarkRed,'LineWidth',1);
hold on
xline(15,'color',colorRfcAwake,'LineWidth',2);
x2 = xline(55,'color',colorRfcNREM,'LineWidth',2);
x1 = xline(97,'color',colorRfcAwake,'LineWidth',2);
xline(105,'color',colorRfcNREM,'LineWidth',2);
x3 = xline(156,'color',colorRfcREM,'LineWidth',2);
xline(224,'color',colorRfcAwake,'LineWidth',2);
xline(248,'color',colorRfcNREM,'LineWidth',2);
xline(342,'color',colorRfcAwake,'LineWidth',2);
xline(360,'color',colorRfcNREM,'LineWidth',2);
xline(450,'color',colorRfcREM,'LineWidth',2);
xline(537,'color',colorRfcAwake,'LineWidth',2);
ylabel('\DeltaD/D (%)')
legend([p3,x1,x2,x3],'Arteriole diameter','Awake','NREM','REM')
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([15,75,135,195,255,315,375,435,495,555,615])
ax34.TickLength = [0.01,0.01];
xlim([15,615])
ylim([-30,60])
% cortical LFP
ax5 = subplot(6,1,5);
semilog_imagesc(T,F,cortNormS,'y')
axis xy
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'Cort LFP','Freq (Hz)'})
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([15,75,135,195,255,315,375,435,495,555,615])
ax5.TickLength = [0.01,0.01];
xlim([15,615])
% hippocampal LFP
ax6 = subplot(6,1,6);
semilog_imagesc(T,F,hipNormS,'y')
axis xy
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (min)')
ylabel({'Hipp LFP','Freq (Hz)'})
set(gca,'box','off')
axis tight
xticks([15,75,135,195,255,315,375,435,495,555,615])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
ax6.TickLength = [0.01,0.01];
xlim([15,615])
% axes properties
ax1Pos = get(ax1,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax5Pos(3:4) = ax1Pos(3:4);
ax6Pos(3:4) = ax1Pos(3:4);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);

function [] = semilog_imagesc(x,y,C,logaxis)
surface(x,y,zeros(size(C)),(C),'LineStyle','none');%
q = gca;
q.Layer = 'top';
if strcmp(logaxis,'y') == 1
    set(gca,'YScale','log');
elseif strcmp(logaxis,'x') == 1
    set(gca,'XScale','log');
elseif strcmp(logaxis,'xy') == 1
    set(gca,'XScale','log');
    set(gca,'YScale','log');
end
axis xy
axis tight
end