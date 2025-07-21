classdef myFigure < handle
    % plot tools for deduction task
    % Zhibing Xiao, Autumn 2023 
    %

    properties 
        DATA
        GROUP
        Colors
        CM
        Width  
        MarkerSize
        FontName
        FontSize
        X
        P
        H
        E
        gca
        g
        enlargeScale
        cmap
    end
 
    properties (Dependent)
        Range
        
    end

    methods

        function obj = myFigure(DATA,GROUP,H, enlargeScale)
            addTools(obj)
            if ~exist('DATA','var'),DATA = [];end
            if ~exist('enlargeScale','var'),enlargeScale=1;end
            obj.enlargeScale =enlargeScale;
          
            if ~exist('GROUP','var')||isempty(GROUP)
                GROUP= ones(size(DATA,1),1);
            end
             
            if ~exist('H','var')||isempty(H),H = obj.setpaper;else,H = obj.setpaper(H);end
            if ~exist('X','var'),X = 1:size(DATA,2);end
            
            obj=obj.set_DATA(DATA);
%             obj.DATA = DATA ;
            obj.H = H;
            obj.GROUP= GROUP;          
            obj.CM = f_deduction_color(1);
            if size(obj.Colors,1)==1,obj.Colors=[.5 .5 .5];end
            
            obj.MarkerSize = 4;
            obj.FontName = 'Arial';
            obj.FontSize = 7;
 
            obj.g = [];
            if isempty(obj.Colors),obj.Colors=obj.CM.gray;end
            obj = f_enlarge(obj);
            obj.cmap = brewermap(64,'Reds');

        end
        function obj = f_enlarge(obj)
            obj.H.Position(3:4) = obj.H.Position(3:4).*obj.enlargeScale;
        end
 
        function [varargout]= addBar(obj,varargin)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            varargin = parseCheckDefault(obj,varargin,'barWidth','Width');
            varargin = parseCheckDefault(obj,varargin,'P','P');
            varargin = parseCheckDefault(obj,varargin,'E','E');
            varargin = parseCheckValue(obj,varargin,'PStarShowNS',false);
            varargin = parseCheckValue(obj,varargin,'PStarBackgroundColor','none');
            varargin = parseCheckValue(obj,varargin,'PLineWidth' ,1);
            varargin = parseCheckValue(obj,varargin,'PLineBackingWidth' ,0.01);
            varargin = parseCheckValue(obj,varargin,'PStarFontSize' ,8);
            varargin = parseCheckValue(obj,varargin,'ErrorbarLineWidth' ,1);
            varargin = parseCheckValue(obj,varargin,'BarLineWidth' ,1);
            
            if ~isempty(obj.E)
%                 obj.E = obj.E(:)';
                offset = max(abs(nanmean(obj.DATA,1))+abs(obj.E))*0.15;
                offset = offset(:)';
                varargin = parseCheckValue(obj,varargin,'PStarOffset',offset );
            end
            if any(mean(obj.DATA,1)<0),estyle = 'I';else,estyle='T';end
            varargin = parseCheckValue(obj,varargin,'ErrorbarStyle' ,estyle);
            
            
            if size(obj.DATA,1)>1
                varargin = parseCheckValue(obj,varargin,'E' ,obj.f_sem(obj.DATA));
            end
   
            m  = nanmean(obj.DATA,1);
 
            figure(obj.H)
            hold on
            h = superbar(obj.X,m,'BarFaceColor',obj.Colors,'PStarShowGT',false,'PStarThreshold',[0.05, 0.01, 0.001], varargin{:});
            xlm = get(gca,'XLim'); %#ok<*CPROPLC> 
            plot(xlm,[0 0],'color','k');%,'LineWidth',0.5);
            xticks(obj.X)
            set(handle(h),'Tag','Bar');
            set(gca,'FontName',obj.FontName);
            set(gca,'FontSize',obj.FontSize);
            obj.gca = gca();            
            varargout{1} =h;
            xlim([obj.X(1)-1 obj.X(end)+1])
             
           
        end
        function varargout=addSpread(obj,varargin)
 
            figure(obj.H)
             
            hold on
            h =plotSpread(obj.DATA,'xValues',obj.X,'distributionColors',obj.Colors*0.7,varargin{:});
            
            set(handle(h{1}),'Tag','Spread');
            varargout{1} =h;
            set(handle(h{1}),'MarkerSize',obj.MarkerSize);
            set(gca,'FontName',obj.FontName);
            set(gca,'FontSize',obj.FontSize);
            
            obj.gca = gca();   
         end
        function varargout=addBarWithSpread(obj,varargin)
            figure(obj.H)                                    
            tmp = addBar(obj,varargin{:});
            width = get(handle(tmp),'BarWidth'); 
            if iscell(width),width=width{1};end

            varargout{1} = addSpread(obj,'xValues',obj.X, ...
                'distributionColors',obj.Colors*0.7,'spreadWidth',width*0.7);
            hold on
            varargout{2} = addBar(obj,varargin{:});
            ylim(obj.Range)
            
            obj = updatePstartLoc(obj);
            obj.gca = gca;
        end
        function varargout=addSpreadWithBar(obj,varargin)
            figure(obj.H)                                    
            tmp = addBar(obj,varargin{:});
            width = get(handle(tmp),'BarWidth');
            if iscell(width),width=width{1};end

%             delete(tmp)

            varargout{1} = addSpread(obj,'xValues',obj.X, ...
                'distributionColors',obj.Colors*0.7,'spreadWidth',width);
            hold on
            varargout{2} = tmp;% addBar(obj,varargin{:});
            ylim(obj.Range)
            
            obj = updatePstartLoc(obj);
            obj.gca = gca;
        end        
        function h = addPstar(obj,msg,P)
            xr = get(gca,'XLim');
            x = mean(xr);
            yr = get(gca,'YLim');
            y  = yr(1)+range(yr)*0.1;
            h=text(x,y,0,msg,"FontSize",obj.FontSize,"FontName",obj.FontName);
            x= (h.Extent(1)+h.Extent(3)*1.1);
            y= (h.Extent(2)+h.Extent(4)*0.6);
            h=text(x,y,0,[f_pValue2flag(P)],'Interpreter','latex','FontSize',obj.FontSize+2);
 
        end
        function h=addScatter(obj,slop,int,varargin)
            if size(obj.Colors,1)~=size(obj.DATA,1)&&size(obj.Colors,1)>0
                C = obj.Colors(1,:);
            else
                C = obj.Colors;
            end
            if isempty(C);C=obj.CM.gray;end
            scatter(obj.DATA(:,1),obj.DATA(:,2),14,C,'filled','o',varargin{:})
%             h=lsline;
            h=[];
            if ~exist('int','var')
                h=lsline;
                h.LineWidth=2;
                h.XData =[min(obj.DATA(:,1)),max(obj.DATA(:,1))];
            elseif ~isempty(slop)
                h=f_plotFittedLine(obj,obj.DATA(:,1),slop,int);
                obj.gca = gca;
            end
            if exist('XLIM','var'),xlim(XLIM);
            else
                xr = obj.DATA(:,1); yr = obj.DATA(:,2);
                XLIM = [min(xr)-0.05*range(xr),max(xr)+0.05*range(xr)];
                YLIM = [min(yr)-0.05*range(yr),max(yr)+0.05*range(yr)];                
                xlim(XLIM);ylim(YLIM);
            end

        end
        function h=addShadeErrorForGLM(obj,weights,xrange,x,y)
            if ~exist('x','var'),x = obj.DATA(:,1);end            
            if ~exist('y','var'),y = obj.DATA(:,2);end            
            if ~exist('weights','var'),weights=ones(size(x));end            
            if ~exist('xrange','var'),xrange=x;end            
            
            md = fitglm(x,y,'Weights',weights);
            newx=linspace(min(xrange),max(xrange),50)';
            [newy,yci]=predict(md,newx,'Alpha',0.05);
            p = polyshape([newx;flip(newx)],[yci(:,1);flip(yci(:,2))]);
    
            hold on
            h=plot(p,'FaceColor',obj.Colors,'EdgeColor','none');            
            alpha(h,0.2)
            plot(newx,newy,'Color',obj.Colors,'LineWidth',1.5)
  
        end
        
 
        function [hs,hl]=addPartialRegressionScatter(obj,XLIM,varargin)
            DATA = obj.DATA;
            DATA(sum(isnan(obj.DATA),2)>0,:)=[];
            x = DATA(:,1);
            y = DATA(:,2);
            c = DATA(:,3:end);
%             [b,ci,stat1]=glmfit(x,c);
%             [b,ci,stat2]=glmfit(y,c);

            [b,ci,stat1]=glmfit(c,x);
            [b,ci,stat2]=glmfit(c,y);            
                        
            xr=stat1.resid+0*median(x);
            yr=stat2.resid+0*median(y);
            corr(xr,yr)
%             hs=scatter(xr,yr,varargin{:});
            col = obj.Colors;
            if size(col,1)==1,
                col = repmat(col,numel(xr),1);
            end
            hs=scatter(xr,yr,15,col ,'filled','CData',col ,varargin{:});
            b = glmfit(xr,yr);r =corr(xr,yr);
            fprintf('r by glm=%.3f, b=%.3f\n',r,b(2));
            if exist('XLIM','var'),xlim(XLIM);
            else
                XLIM = [min(xr)-0.05*range(xr),max(xr)+0.05*range(xr)];
                YLIM = [min(yr)-0.05*range(yr),max(yr)+0.05*range(yr)];                
                xlim(XLIM);ylim(YLIM);
            end
%             hl=lsline;
%             hl.LineWidth=2;
            rng= [min(xr),max(xr)];
            hold on
            obj.addShadeErrorForGLM(ones(size(xr)),rng,xr,yr)

            obj.gca = gca;
        end        
                
        function h =addShadeErrorBar(obj,x,y,e,varargin)
            figure(obj.H)
            if ~exist('x','var'),x=obj.X;end
            if ~exist('y','var'),y=nanmean(obj.DATA,1);end
            if ~exist('e','var'),e=obj.f_sem(obj.f_sem(obj.DATA));end
%             c = obj.alphaColor(obj.Colors,0.2);
            % varargin = parseCheckValue(obj,varargin,'color',obj.Colors);

            if ~isempty(varargin) 
                if numel(varargin)==1
                    lineprop = varargin{1};
                    trans=[];
                else
                    lineprop = varargin;
                    trans =varargin{end};
                end
            else
                lineprop = {};
                trans =[];
            end

            h=shadedErrorBar(x,y,e,lineprop,trans);
            h.mainLine.LineWidth
%             h.mainLine.LineWidth= 2;
            obj.gca = gca;
        end
        function se = f_sem(obj,d)
            if nargin<2,d = obj.DATA;end
            if size(d,1)>1
                se = nanstd(d,1)./(sqrt(sum(~isnan(d),1)));
            else
                se = [];
            end
%             if isempty(obj.E),obj.E = se;end
        end

        function obj = updatePstartLoc(obj,force,rate)
            if ~exist('rate','var'),rate=0.98;end
            if ~exist('force','var'),force=false;end
            if min(size(obj.P))>1&&~force,
                CASE = 1;
            else
                CASE = 2;
                return;
            end
            hc=get(obj.gca,'Children');
            
            ht=hc(ismember(get(hc,'Type'),'text'));
            
            if ~isempty(ht) && CASE==1
                hs=ht(containstr(obj,get(ht,'String'),{'∗' 'n.s.'}));
                if ~isempty(hs)
                    ps =cat(1,hs(:).Position);
                    if any(ps(:,2)>obj.gca.YLim(2))
                        dev = max(ps(:,2)-obj.gca.YLim(2));
                        for i = 1: numel(hs)
                            hs(i).Position(2) =hs(i).Position(2)-dev;
                        end
                    end
                end
            end
            if ~isempty(ht) && CASE==2
                hs=ht(containstr(obj,get(ht,'String'),{'∗' 'n.s.'}));
                for i = 1: numel(hs)
                    hs(i).Position(2) =obj.gca.YLim(2)*rate;
                end
            end            
           
            obj.gca = gca();
        end
        
        function ylim(obj,limit)
            if ~exist('limit','var')
                set(obj.gca,"YLim",obj.Range*0.95)
            elseif numel(limit)>1
                set(obj.gca,"YLim",limit)
            else
                set(obj.gca,"YLim",obj.Range*limit)
            end
            obj = updatePstartLoc(obj);
        end        
        function color = alphaColor(~,color,alpha)
            c = rgb2hsv(color);
            c(2) = c(2)*(alpha);
            c(3) = c(3)*(1-alpha);
            color= hsv2rgb(c);
        end
        function savefig(obj,sfilename,H)
            if nargin<3,H = obj.H;end
            fprintf('saving current figure...')
            if nargin<1,sfilename = fullfile(pwd,'figure.fig');end
            [d,f,~] = fileparts(sfilename);
            sfilename = fullfile(d,[f,'.pdf']);
            if ~exist(d,'dir')&&~isempty(d)
                mkdir(d);end

            print(H,sfilename,'-dpdf', '-vector')
%             savefig(H,strrep(sfilename,'.pdf','.fig'))
            exportgraphics(H,strrep(sfilename,'.pdf','.png'),'Resolution',600)
            fprintf('saved\n')

        end

        function [h_leg,h_leg_icon]=legend(obj,handels,LEGENDS,varargin)
            
            varargin = parseCheckValue(obj,varargin,'location','eastoutside');
            [h_leg,h_leg_icon] = legend(handels,LEGENDS, ...
                'fontSize',obj.FontSize,'fontName',obj.FontName, ...
                varargin{:});
            pos = h_leg.Position;
%             hw_rate = pos(4)/pos(3);
%             h_leg.Position = [pos(1)*0.5 0.75 pos(3:4)];
            h = findobj(h_leg_icon,'type','patch');
            width  = range(h(1).XData) ;
            hw_rate = (range(h(1).YData)*pos(4))/(range(h(1).XData)*pos(3)) ;
%             for ii = 1:length(h), if length(h(ii).XData)==2, h(ii).XData(1) = h(ii).XData(1) + 0.5*range(h(ii).XData); end; end
            for ii = 1:length(h) 
                h(ii).XData(3:4)=h(ii).XData(1:2)+1.6*width*hw_rate;
            end
            h = findobj(h_leg_icon,'type','text');
            for ii = 1:length(h)
                h(ii).Position(1)=h(ii).Position(1)-(2.5*width*hw_rate );
            end
            legend boxoff
%             h_leg.Position
            obj.setFont
            disp('saved')
        end

        function setFont(obj,axs)
            if ~exist('axs','var');axs = findall(obj.H,'type','axes');
            end
            for ax = axs                
                set(ax,'FontName',obj.FontName);
                set(ax,'FontSize',obj.FontSize-1);
                set(ax,'Layer','top','Color','none')
                ht=findall(ax,'Type','Text');      
                hs=containstr(obj,get(ht,'String'),{'∗' '*' 'n.s.'});
                 
             
                set(ht(~hs),'FontName',obj.FontName);
                set(ht(~hs),'FontSize',obj.FontSize);           
            end
        end
        function H=setpaper(obj,H)
%                H = figure('Position',[1400 100 800 800],'Color','w');
%               H = figure('Position',[1400 100 1000 1000],'Color','w'); % A4 size
%               H = figure('Position',[1400 100 2000 1000],'Color','w'); % A4 size
             if ~exist('H','var')&& ~isempty(groot) && ~isempty(groot().Children) && isvalid(groot().CurrentFigure)
                H = gcf; 
                if isempty(H.Children) && strcmpi(H.Tag,'cFig') % use old blank
                    figure(H)
                    return
                end        
             end
             
             if exist('H','var')&& ishandle(H)
                 H.Units = "centimeters";
                 if ~isequal(H.Position(3:4),[18 18])
                    H.Position = [50 2 18*obj.enlargeScale 18*obj.enlargeScale];
                 end
                 H.Color=[1 1 1];
                 H.Tag = 'cFig'; 
             else
                 H = figure('unit','centimeters', ...
                     'Position',[50 2 18*obj.enlargeScale 18*obj.enlargeScale],'Color','w','Tag','cFig'); % same for visulization, print,and adoube illustrator
             end
             H.PaperUnits = 'centimeters';
             H.PaperPosition(1:4)=[0 29-18*obj.enlargeScale 18*obj.enlargeScale 18*obj.enlargeScale];
             H.PaperType='A4';
             figure(H)
        end
        function h=f_plotFittedLine(obj,x,slop,int)            
            xx = [min(x),max(x)];
            axes(gca),hold on;
            h  = [];
            for ii =1:numel(slop)
                yy = xx*slop(ii)+int(ii);
                h=plot(xx,yy);
                h.LineWidth=2;
                h.Color =obj.Colors;
                h(ii) = h;
            end
        end
        function h = addHhatchfill2(obj,barh,type,varargin)
            varargin = obj.parseCheckValue(varargin,'HatchAngle',45);
            varargin = obj.parseCheckValue(varargin,'hatchcolor',obj.CM.gray);
            hatchfill2(barh,type,varargin{:})
            obj.gca = gca;
        end
        function h = addPatch(obj,X,Y,Z,C,varargin)
            
        end

        function P = f_symmetricMatrix(obj,P)
                pp = nan([size(P),2]); pp(:,:,1) = P; pp(:,:,2)=P';
                P =squeeze(nanmean(pp,3));
        end

        function ha = addSubplotName(obj,str,r,c)
            border =0.13; width=(1-2*border)/r; height = (1-2*border)/c;
            vx = border:(1-2*border)/r:1-border ;
            vy = flip(border:(1-2*border)/c:1-border) ;
            [~,ix]=min(abs(vx-obj.gca.Position(1)));
            [~,iy]=min(abs(vy-obj.gca.Position(2)));
%             [~,iy]=min(abs(vy-sum(obj.gca.Position([2]))));
%             obj.gca.YLabel.Position(1)
            loc = [vx(ix)-0.1*width, vy(iy)+1.2*height,0.01 0.01];
            ha=annotation('textbox',loc,'String',str,'EdgeColor','none','FontSize',16,'FontName','Arial');
        end
        function ha = addGcaName(obj,str)
            border =0.13; width=obj.gca.XLim; height = obj.gca.YLim;
 
            loc = [width(1)-0.3*range(width), height(2)+0.1*range(height),0.01 0.01];
%             ha=text('textbox',loc,'String',str,'EdgeColor','none','FontSize',16,'FontName','Arial');
             ha=text(loc(1),loc(2),-1,str,'FontSize',16,'FontName','Arial');
        end    
        function  varargout = f_plotSurf(obj, cfg, ax)    
            if ~exist("cfg",'var'),cfg =struct;end
            if ~exist("ax",'var'),ax =obj.gca;end
            cpath = path;
            addSurfPlot(obj,true) % add path 
            
            cfgDefault = defaultSurfCfg(obj);
            cfg.axis = ax; axes(ax); %set(gca,'units','pixels');


            opts_def = fieldnames(cfgDefault);
            opts_inp = fieldnames(cfg);
            for i = 1:numel(opts_inp)
                copt = opts_inp{i};
%                 if any(ismember(opts_def,copt))
%                     if ~isequal(cfgDefault.(copt),cfg.(copt))
%                         cfgDefault.(copt)=cfg.(copt);
%                     end
%                 else
                    cfgDefault.(copt)=cfg.(copt);
%                 end
            end

            cfgOut=plotPialSurf('fsaverage',cfgDefault);
            obj = pullOutElect(obj,ax,cfgOut,cfgDefault.elecColorScale);
            varargout{1} = cfgOut; 
            restoredefaultpath % remove path           
            addpath(cpath)
        end
        function h=f_plotSurfOnHipp(obj,CM,channels,value,clim)
            if ~exist("CM",'var')
                paths.task_lib = fullfile('/HDD/lbcode-dev-xzb/analysis/zhibing/deduction/ripple','local_task_lib');
                CM = load(fullfile(paths.task_lib,'mycmappos.mat'));
                CM = CM.mycmap;
            end
            hd=load('/HDD/Results/manual/hipp_visualization.mat');
            includedRows = find(hd.elecinfo.Helec>0);
%             Helec = hd.elecinfo.Helec(includedRows);
            Helec = hd.Helec ;  h = hd.Hes;
            if ~exist('clim','var')
                clim=[0.0 0.6]; 
            end
            if ~exist('thr ','var')
                thr = 0.0;
            end            
            value(abs(value)<thr) = thr; value(value>clim(2))=clim(2); value(value<clim(1))=clim(1);  
            grouplabel = hd.elecinfo.groupLabels(includedRows);
            [~,ia,ib]=intersect(grouplabel,channels,'stable');
            v = nan(size(grouplabel));
            v(ia) = value(ib);           
            axs=findall(h,'type','axes');
            axs(1).Position([1 2])=[0.02 0.000];
            axs(2).Position([1 2])=[0.4 0.000];
            
             
            ncolors = size(CM,1);
            for ii = 1:size(v,1)
                colind = floor(interp1(clim,[1 ncolors],v(ii),'linear',ncolors));
                set(Helec(ii),'facecolor',CM(colind,:));
            end
            axes('Position',[0.8 0.2 0.02 0.6]);axis off
            cb=colorbar('vertical');
%             cb = cbar('vert')
            caxis(clim); colormap(CM);
            cb.Label.String = 'value';cb.Label.FontSize=8;
            cb.Limits = clim; cb.Ticks = clim;
            cb.Position(3)=0.02;
        end
        function [g, obj]= gramm(obj,varargin)
            x = repmat(1:size(obj.DATA,2),size(obj.DATA,1),1);
            x = x(:);
            y = obj.DATA(:);
            if isempty(obj.Colors)
                color = x;
            else
                c = arrayfun(@(x) {obj.Colors(x,:)},1:size(obj.Colors,1),'un',0);
                c = reshape(c,1,numel(c));
                color = repmat(c,size(obj.DATA,1),1);
                color = cat(1,color{:});
                color = cell2mat(color);
%                 color = repmat(obj.Colors,)
            end

            varargin = parseCheckValue(obj,varargin,'x',x);
            varargin = parseCheckValue(obj,varargin,'y',y);
            varargin = parseCheckValue(obj,varargin,'color',color);
%             g = gramm('x',cars.Origin_Region,'y',cars.Horsepower,'color',cars.Origin_Region,'subset',cars.Cylinders~=3 & cars.Cylinders~=5);
            g = gramm(varargin{:});

            obj.g = g;
        end
        function obj = gramm_setloc(obj,pos)
           
            loc_facet = obj.g.facet_axes_handles.Position;
%             axs=fieldnames(obj.g.axe_property);

            set(obj.g.facet_axes_handles,'Position',pos);

            loc_facet_new = obj.g.facet_axes_handles.Position;
            factor = loc_facet_new-loc_facet;

            loc_l = [loc_facet_new(1)+loc_facet_new(3),loc_facet_new(2),loc_facet_new(3),loc_facet_new(4)];
            obj.g.legend_axe_handle.Position=loc_l;
 
        end
 
        function lb = f_DMNsubNames(obj,sub2query,inshort,forceInshort)
            DMNsubs = {'DMNmPFC'  'DMNPCC'  'DMNTPJ' 'DMNTL' 'DMNLOFC' 'DMNDLPFC'...
                 'Y7_Visual' 'Y7_Default' 'Y7_DorsalAttention' 'Y7_Frontoparietal' ...
                'Y7_VentralAttention' 'Y7_Limbic' 'Y7_Somatomotor'};
            DMNlabels = {'mPFC' 'PCC'  'Inferior parietal' 'Middle temporal' 'Inferior frontal' 'Superior frontal'...
                'Visual' 'Default' 'Dorsal Attention' 'Frontoparietal' 'Ventral Attention' 'Limbic' 'Somatomotor'};
            ShortCut = {'mPFC' 'PCC'  'Inferior parietal' 'Middle temporal' 'Inferior frontal' 'Superior frontal'...
                'VIS' 'DMN' 'DA' 'FPN'  'VA' 'LN' 'SM'};
            ShortCut_force = {'mPFC' 'PCC'  'IP' 'MT' 'IF' 'SF'...
                'VIS' 'DMN' 'DA' 'FPN'  'VA' 'LN' 'SM'};            
            if exist('inshort','var')&&inshort
                DMNlabels = ShortCut;
            end
            if exist('forceInshort','var')&&forceInshort
                DMNlabels = ShortCut_force;
            end
            
            
            [~,ia,ib]=intersect(sub2query,DMNsubs,'stable');
            lb = DMNlabels(ib);
            
        end 
        function obj=set_DATA(obj,x)
            obj.DATA = x;
            obj.E = obj.f_sem(x);
            if size(x,2)~=numel(obj.P)
                obj.P = [];
            end
        end        
        function set.P(obj,x)
%             obj.Range = x ;
            tmp = x; tmp(isnan(tmp))=0;
            if all(size(x)>1)&& ~issymmetric(tmp)
                x = obj.f_symmetricMatrix(x);
            end
            obj.P = x;
        end
        function set.H(obj,x)
            obj.H = obj.setpaper(x);
        end
        
        function x = get.Range(obj)
            x =  [min(prctile(obj.DATA,5)), max(prctile(obj.DATA,95))];
        end
        function x = get.E(obj)
            x =  obj.f_sem(obj.DATA);
        end
        function set.DATA(obj,x)
%             if size(x,1)>1&&size(x,2)==1,x = x';end
            obj.DATA = x;
            updateProperties(obj);
        end

        function obj=updateProperties(obj)
            dimDATA = size(obj.DATA,2);
            if numel(obj.X)~=dimDATA 
                obj.X = 1:size(obj.DATA,2);
            end      
            if numel(obj.P)~=dimDATA
                obj.P = [];
            end    
            if size(obj.Colors,1)~=dimDATA 
                obj.Colors = brewermap(dimDATA,'RdYlGn');
            end
            

        end
        function pulledFigureCrd=pullElectrodesTowardsTheCamera(obj,figureCrd,axisHandle,pullingExtent)

            cameraVec=get(axisHandle,'CameraPosition')-get(axisHandle,'CameraTarget');
            normedCameraVec=cameraVec/norm(cameraVec);
            pulledFigureCrd=bsxfun(@plus,figureCrd,normedCameraVec*pullingExtent);
        end
        
     end


    %%
    methods (Access=private)
%         function ind = containstr(scell,str)
%             ind = false(sizes(scell));
%             for ii = 1:numel(scell)
%                 s = scell{ii};
%                 if ~isempty(s)
%                     if ischar(s)||isstring(s)
%                         ind(ii)=contains(s,str);
%                     elseif iscell(s)
%                         ind(ii)=any(contains(s,str));
%                     end
%                 end
%             end
%         end
        function obj = pullOutElect(obj,ax,cfgOut,scale)
            axes(ax);
            allElec = findobj(gca,'tag','Electrode');
%             set(allElec,'LineWidth',0.5,'markeredgecolor',COLOR.gray);
%             GrayElec = findobj(gca,'markerfacecolor',COLOR.belowThr);
            nonGrayElec = allElec;%setdiff(allElec,findobj(gca,'markerfacecolor',COLOR.belowThr),'stable');
            cols =cat(1,nonGrayElec.MarkerFaceColor);
            [c,ord]=sort(sum(cols,2),1,'ascend');
            nonGrayElec = nonGrayElec(ord);
            if ~isempty(nonGrayElec)
                X = [nonGrayElec.XData];
                Y = [nonGrayElec.YData];
                Z = [nonGrayElec.ZData];
                S = [nonGrayElec.MarkerSize];
                pullingExtent = 5;
                maxEnlargingExtent = 1;
                cmap = cfgOut.elecCmapName;
                ncolors = size(cmap,1)-1;
                pulledCrd = [];
                for i = 1:length(nonGrayElec)
                    current_color = [nonGrayElec(i).MarkerFaceColor];
                    ii = find(mean(bsxfun(@minus,cmap,current_color),2)==0)-1;
                    thrdist =0;% (scale(3)/scale(2))*0.5;
                    enlarging_factor = (abs((ii/ncolors)-0.5)-thrdist)/0.5; % range: [0, 1];
                    pulledCrd=pullElectrodesTowardsTheCamera([X(i),Y(i),Z(i)],gca,pullingExtent);
                    set(nonGrayElec(i),'XData',pulledCrd(1),'YData',pulledCrd(2),'ZData',pulledCrd(3));%,'markersize',S(i)*(1+enlarging_factor*maxEnlargingExtent));
%                     set(nonGrayElec(i),'markeredgecolor','k','LineWidth',0.5)
                end
            end
        end
        function obj=setWidth(obj)
            hc=get(obj.gca,'Children');            
        end
        function ind=containstr(obj,str,s)
            ind=false(numel(str),1);
            if isempty(str),return;end
            if ~iscell(str)
                ind=contains(str,s);
                return
            end
            for i = 1:numel(str)      
                if ischar(str{i})
                    if any(contains(str(i),s))&&~isempty(str(i))
                        ind(i) = true;
                    end
                end
            end
        end
        function addTools(obj)
            addpath('/HDD/ripple-representation-learning/toolboxs/brewermap')
            addpath('/HDD/ripple-representation-learning/toolboxs/export_fig-3.34')
            addpath('/HDD/ripple-representation-learning/toolboxs/plotSpread')
            addpath('/HDD/ripple-representation-learning/toolboxs/scottclowe-superbar-ce333a9/superbar')
  %        addpath('/HDD/lbcode-dev-xzb/analysis/zhibing/deduction/utils/toolbox/hatchfill2/')
    %         addpath('/HDD/lbcode-dev-xzb/analysis/zhibing/deduction/utils/toolbox/gramm-master')
        end
        function addSurfPlot(obj,add)
            rootdir = '/HDD';
            if add
                addpath(fullfile(rootdir,'lbcode-dev-xzb/analysis/zhibing/deduction/ripple/local_lib'))
                addpath(genpath(fullfile(rootdir,'lbcode-dev-xzb/externalPackages/IELVis/')))
                addpath(genpath(fullfile(rootdir,'lbcode-dev-xzb/externalPackages/afni_matlab/')))
                addpath(genpath(fullfile(rootdir,'lbcode-dev-xzb/externalPackages/freezeColors/')))
                addpath(fullfile(rootdir,'lbcode-dev-xzb/externalPackages/brewermap/'))
            else
                rmpath(fullfile(rootdir,'lbcode-dev-xzb/analysis/zhibing/deduction/ripple/local_lib'))
                rmpath(genpath(fullfile(rootdir,'lbcode-dev-xzb/externalPackages/IELVis/')))
                rmpath(genpath(fullfile(rootdir,'lbcode-dev-xzb/externalPackages/afni_matlab/')))
                rmpath(genpath(fullfile(rootdir,'lbcode-dev-xzb/externalPackages/freezeColors/')))
                rmpath(fullfile(rootdir,'lbcode-dev-xzb/externalPackages/brewermap/'))
            end
        end
        function [inputVars,v]=parseCheckDefault(obj,inputVars,target,propertTag)
            % check if the 'target' in inputVars match the obj.(propertTag)
            % if not, overwirte obj.(propertTag) by target
            % if target not exist, and obj.(propertTag) is not empty, add
            % it to pararms 
            if ~isempty(inputVars)
                params = [inputVars(1:2:end)',inputVars(2:2:end)',];
                ind=ismember(params(:,1),target);
                if sum(ind)>1,error('duplicated pararmeters');end
                if any(ind)                
                    attrib = params{ind,2};
    %       
                        if ~isequal(attrib,obj.(propertTag))
                            obj.(propertTag) = attrib;
                        end
                else
                    if ~isempty(obj.(propertTag))
                        inputVars=[inputVars,{target},{obj.(propertTag)}];
                    end
                end
            else
                if ~isempty(obj.(propertTag))
                    inputVars=[inputVars,{target},{obj.(propertTag)}];
                end
            end
            v = obj.(propertTag);
 

        end     
        function inputVars=parseCheckValue(obj,inputVars,target,tarValue)
            % check if the 'target' in inputVars match the obj.(propertTag)
            % if not, overwirte obj.(propertTag) by target
            % if target not exist, and obj.(propertTag) is not empty, add
            % it to pararms 
            if ~isempty(inputVars)
                params = [inputVars(1:2:end)',inputVars(2:2:end)',];
                ind=ismember(lower(params(:,1)),lower(target));
                if sum(ind)>1,error('duplicated pararmeters');end
                if any(ind)                
                    attrib = params{ind,2};
                    if ~isequal(attrib,tarValue)
                        tarValue = attrib;
                    end
                end
            end
            fprintf('param %s added:\n',target)
            disp(tarValue)
            inputVars=[inputVars,{target},{tarValue}];
        end          
        function cfg=defaultSurfCfg(obj)        
            cfg = struct;
            cfg.figId = obj.H;
            cfg.clearFig = 'n';
            cfg.elecCbar = 'y';
            % cfg.elecColorScale= scale;
            cfg.showLabels='n';
            cfg.surfType = 'inflated';
            cfg.edgeBlack = 'n';
            cfg.elecShape = 'marker';
            % cfg.elecShape ='sphere';
            % cfg.elecUnits = cbunits;
            cfg.elecSize = 3;
            cfg.opaqueness = 1;
            cfg.pullOut = 5;
            cfg.overlayParcellation = 'Y7';
            cfg.title='';
            cfg.elecCmapName = 'mycmap';
            cfg.overlayParcellation='Y7';
            cfg.parcellationColors = ones(8,3)*187;%  *0.8*255;
            cfg.parcellationColors(8,:)=[190 161 165]; %[244 196 197] % [0.8 0.5 0.1]*256
            cfg.elecCbar = 'y';
            cfg.view = 'lm';

            
            
        end

    end

end
