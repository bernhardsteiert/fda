close all; clear all;figure('Position',[50 50 350 800]);subplot1(4,2,'Min',[0.01 0.01],'Max',[0.99 0.99],'Gap',[0.01 0.01]);
subplot1(1);
EGF_parental = imread('.\Workspaces\parental-15MIN\007007-9-001001003.TIF','PixelRegion',{[117 117+615-1],[385 385+615-1]});
imshow(imadjust(EGF_parental,[0 0.013],[0.03 1])); title('Parental: EGF 100 ng/ml'); axis image;

subplot1(2);
EGF_reporter = imread('.\Workspaces\reporter-15MIN\007007-5-001001002.TIF','PixelRegion',{[1 615],[25 25+615-1]});
imshow(imadjust(EGF_reporter,[0 0.01],[0.05 1])); title('Reporter: EGF 100 ng/ml'); axis image;


subplot1(3);
EGF_parental = imread('.\Workspaces\parental-15MIN\007009-12-001001003.TIF','PixelRegion',{[117 117+615-1],[385 385+615-1]});
imshow(imadjust(EGF_parental,[0 0.013],[0.03 1])); title('Parental: EGF 2 ng/ml'); axis image;

subplot1(4);
EGF_reporter = imread('.\Workspaces\reporter-15MIN\007009-10-001001002.TIF','PixelRegion',{[1 615],[25 25+615-1]});
imshow(imadjust(EGF_reporter,[0 0.01],[0.05 1])); title('Reporter: EGF 2 ng/ml'); axis image;


subplot1(5);
EGF_parental = imread('.\Workspaces\parental-15MIN\007011-6-001001003.TIF','PixelRegion',{[150 150+615-1],[170 170+615-1]});
imshow(imadjust(EGF_parental,[0 0.025],[0.1 1])); title('Non-stimulated'); axis image;

subplot1(6);
EGF_reporter = imread('.\Workspaces\reporter-60MIN\002011-1-001001002.TIF','PixelRegion',{[350 350+615-1],[1 1+615-1]});
imshow(imadjust(EGF_reporter,[0 0.015],[0.05 1])); title('Non-stimulated'); axis image;

subplot1(7);
EGF_parental = imread('.\Workspaces\parental-15MIN\002011-6-001001003.TIF','PixelRegion',{[150 150+615-1],[170 170+615-1]});
imshow(imadjust(EGF_parental,[0 0.025],[0.05 1])); title('MK2206 1\muM'); axis image;

subplot1(8);
EGF_reporter = imread('.\Workspaces\reporter-15MIN\002011-10-001001002.TIF','PixelRegion',{[150 150+615-1],[580 580+615-1]});
imshow(imadjust(EGF_reporter,[0 0.02],[0.05 1])); title('MK2206 1\muM'); axis image;
