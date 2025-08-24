classdef myFigure < handle
    % plot tools  
    % Author, Zhibing Xiao, Autumn 2023 
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
        VarName
        
    end
 
    properties (Dependent)
        Range
        
    end
    properties(Access=private)
        HsgTitle 
    end
    methods

        function obj = myFigure(DATA,GROUP,H)
            addTools(obj)
            if ~exist('DATA','var'),DATA = [];end
          
            if ~exist('GROUP','var')||isempty(GROUP)
                GROUP= ones(size(DATA,1),1);
            end
             
            if ~exist('H','var')||isempty(H),H = obj.setpaper;else,H = obj.setpaper(H);end
            if ~exist('X','var'),X = 1:size(DATA,2);end
            
            obj.DATA = DATA;
            obj.H = H;
            obj.GROUP= GROUP;          
            obj.CM = f_deduction_color(1);
            if size(obj.Colors,1)==1,obj.Colors=[.5 .5 .5];end
            
            obj.MarkerSize = 4;
            obj.FontName = 'Arial';
            obj.FontSize = 7;
 
            obj.g = [];
            if isempty(obj.Colors),obj.Colors=obj.CM.gray;end
            obj.cmap = brewermap(64,'Reds');

        end
 
        function addLineForPairedSpread(obj,varargin)
            xx = repmat(obj.X,height(obj.DATA),1);
            yy = obj.DATA;
            plot(xx',yy',varargin{:})
        end

        function x = addSgtitle(obj,text)
            obj.H.Name = text;
            delete(handle(obj.HsgTitle))
            text = strrep(text,'_','\_');
            x=annotation(obj.H,'textbox',[0 0.98 1 0.02],'String',text, ...
                'HorizontalAlignment','center','EdgeColor','none','FontWeight','bold','FontName',obj.FontName,'FontSize',obj.FontSize+2,'Tag','sgtitle');
            obj.HsgTitle = x;
            obj.setFont();
        end        
        function [varargout]= addBar(obj,varargin)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            % set some defaults
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
            varargin = parseCheckValue(obj,varargin,'BarEdgeColor' ,[1 1 1]*0.7);            
            if ~isempty(obj.E)
%                 obj.E = obj.E(:)';
                offset = max(abs(nanmean(obj.DATA,1))+abs(obj.E))*0.15;
                offset = offset(:)';
                varargin = parseCheckValue(obj,varargin,'PStarOffset',offset );
            end
            if any(mean(obj.DATA,1)<0),estyle = 'I'; else, estyle='T';end
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
            x = max(xr);
            yr = get(gca,'YLim');
            y  = yr(1)+range(yr)*0.1;

            if contains(msg,'r')
                 msg = strrep(msg,msg(1),['\it{' msg(1) '} \rm ']);
                 msg = strrep(msg,'  ',' ');
            end
            pflag = f_pValue2flag(P);
            if P>0.05,pflagsize = obj.FontSize+1;else,pflagsize = obj.FontSize+3;end
            h=text(x,y,0,pflag,'FontSize',pflagsize,'HorizontalAlignment','right','VerticalAlignment','middle');
            h=text(x-h.Extent(3),y,0,[msg ],'FontSize',obj.FontSize ,'HorizontalAlignment','right');
     
        end
        function h=addScatter(obj,slop,int,varargin)
            if size(obj.Colors,1)~=size(obj.DATA,1)&&size(obj.Colors,1)>0
                C = obj.Colors(1,:);
            else
                C = obj.Colors;
            end
            if isempty(C);C=obj.CM.gray;end
            scatter(obj.DATA(:,1),obj.DATA(:,2),35,C,'filled','o',varargin{:},'MarkerEdgeColor',[1 1 1])
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
      
            xr = obj.DATA(:,1); yr = obj.DATA(:,2);
            XLIM = [min(xr)-0.05*range(xr),max(xr)+0.05*range(xr)];
            YLIM = [min(yr)-0.05*range(yr),max(yr)+0.05*range(yr)];
            xlim(XLIM);ylim(YLIM);


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
            hs=scatter(xr,yr,30,col ,'filled','CData',col ,varargin{:},'MarkerEdgeColor',[1 1 1]);
            b = glmfit(xr,yr);r =corr(xr,yr);
            fprintf('r by glm=%.3f, b=%.3f\n',r,b(2));
            if exist('XLIM','var'),xlim(XLIM);
            else
                XLIM = [min(xr)-0.05*range(xr),max(xr)+0.05*range(xr)];
                YLIM = [min(yr)-0.05*range(yr),max(yr)+0.05*range(yr)];                
                xlim(XLIM);ylim(YLIM);
            end
 
            rng= [min(xr),max(xr)];
            hold on
            obj.addShadeErrorForGLM(ones(size(xr)),rng,xr,yr)

            obj.gca = gca;
        end        
                
        function h =addShadeErrorBar(obj,x,y,e,varargin)
            figure(obj.H)
            if ~exist('x','var'),x=obj.X;end
            if ~exist('y','var'),y=nanmean(obj.DATA,1);end
            if ~exist('e','var'),e=obj.f_sem(obj.DATA);end
 
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
            if nargin<2,sfilename = fullfile(pwd,'figure.fig');end
            [d,f,~] = fileparts(sfilename);
            sfilename = fullfile(d,[f,'.pdf']);
            if ~exist(d,'dir')&&~isempty(d)
                mkdir(d);end

            print(H,sfilename,'-dpdf', '-vector')
            exportgraphics(H,strrep(sfilename,'.pdf','.png'),'Resolution',300)
            fprintf('saved\n')

        end

        function [h_leg,h_leg_icon]=legend(obj,handels,LEGENDS,varargin)
            
            varargin = parseCheckValue(obj,varargin,'location','eastoutside');
            [h_leg,h_leg_icon] = legend(handels,LEGENDS, ...
                'fontSize',obj.FontSize,'fontName',obj.FontName, ...
                varargin{:});
            pos = h_leg.Position;
            h_leg.ItemTokenSize = [30 9];
            
            legend boxoff
 
            obj.setFont
            disp('saved')
        end

        function setFont(obj,axs)
            if ~exist('axs','var')||isempty(axs);axs = findall(obj.H,'type','axes');
            end
            for ax = axs                
                set(ax,'FontName',obj.FontName);
                set(ax,'FontSize',obj.FontSize);
                set(ax,'Layer','top','Color','none')
                ht=findall(ax,'Type','Text');      
                hs=containstr(obj,get(ht,'String'),{'∗' '*' 'n.s.'});
                 
             
                set(ht(~hs),'FontName',obj.FontName);
                set(ht(~hs),'FontSize',obj.FontSize);           
            end
        end
        function H=setpaper(obj,H)
 
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
                    H.Position = [50 2 18 18];
                 end
                 H.Color=[1 1 1];
                 H.Tag = 'cFig'; 
             else
                 H = figure('unit','centimeters', ...
                     'Position',[50 2 18 18],'Color','w','Tag','cFig'); % same for visulization, print,and adoube illustrator
             end
             H.PaperUnits = 'centimeters';
             H.PaperPosition(1:4)=[0 29-18 18 18];
             H.PaperType='A4';
             figure(H)
        end
        function h=f_plotFittedLine(obj,x,slop,int)            
            xx = [min(x),max(x)];
            axes(gca),hold on;
            hh  = [];
            for ii =1:numel(slop)
                yy = xx*slop(ii)+int(ii);
                h=plot(xx,yy);
                h.LineWidth=2;
                h.Color =obj.Colors;
                hh(ii) = h;
            end
        end
 

        function P = f_symmetricMatrix(obj,P)
                pp = nan([size(P),2]); pp(:,:,1) = P; pp(:,:,2)=P';
                P  = squeeze(nanmean(pp,3));
        end

        function ha = addSubplotName(obj,str,r,c)
            border =0.13; width=(1-2*border)/r; height = (1-2*border)/c;
            vx = border:(1-2*border)/r:1-border ;
            vy = flip(border:(1-2*border)/c:1-border) ;
            [~,ix]=min(abs(vx-obj.gca.Position(1)));
            [~,iy]=min(abs(vy-obj.gca.Position(2)));
            loc = [vx(ix)-0.1*width, vy(iy)+1.2*height,0.01 0.01];
            ha=annotation('textbox',loc,'String',str,'EdgeColor','none','FontSize',16,'FontName','Arial');
        end
        function ha = addGcaName(obj,str)
            border =0.13; width=obj.gca.XLim; height = obj.gca.YLim;
 
            loc = [width(1)-0.3*range(width), height(2)+0.1*range(height),0.01 0.01];
%             ha=text('textbox',loc,'String',str,'EdgeColor','none','FontSize',16,'FontName','Arial');
             ha=text(loc(1),loc(2),-1,str,'FontSize',16,'FontName','Arial');
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
    
        function set.DATA(obj,x)
%             if size(x,1)>1&&size(x,2)==1,x = x';end
            if istable(x)
                obj.VarName = x.Properties.VariableNames;
                x = x.Variables;
            end
            obj.DATA = x;
            obj.E = obj.f_sem(x);
            if size(x,2)~=numel(obj.P)
                obj.P = [];
            end
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
    
        % 
     end

methods (Static)

           function [hs,rg]=shadeErrorBar_matlab(x,data,edgeOrError,color,ax)
            
                if nargin==1
                    data = x;
                    x = 1:size(data,2);
                end
                if ~exist('color','var'),color = [0.3 0.3 0.3]; end
                if ~exist('ax','var'),ax = gca; end
                data = data(~all(isnan(data),2),:);
                m = mean(data,1);
                if ~exist('edgeOrError','var')||isempty(edgeOrError), edgeOrError = std(data,1)./sqrt(size(data,1)); end
                if width(edgeOrError)==1||height(width(edgeOrError))
                    edgeOrError = edgeOrError(:);
                    edgeOrError = [m'-edgeOrError, m'+edgeOrError];
                end
                idx_nan =isnan(m);
                if size(edgeOrError,1)<size(edgeOrError,2), edgeOrError =edgeOrError';end
                wx = x;
                hs = [];
                
                sx = [wx(~idx_nan ) flip(wx(~idx_nan )) ];
                sy = [edgeOrError(~idx_nan ,1) ;flip(edgeOrError(~idx_nan ,2))];
                hs.shade = fill(ax,sx,sy, color,'FaceAlpha',0.1,'EdgeColor',[1 1 1]*0.9); hold on
                hs.mainLine = plot(ax,wx,m,'-k','LineWidth',1.5,'Color',color);hold on                
                rg = edgeOrError(:);
                rg = [min(rg(~isnan(rg))), max(rg(~isnan(rg)))];
                ylim(rg*1.1)
                xlim([x(1) x(end)]);
                

        end
         function subplot_compress_margins(left_margin, right_margin)
            if nargin < 1
                left_margin = 0.05;
            end
            if nargin < 2
                right_margin = 0.05;
            end

            all_axes = findall(gcf, 'Type', 'axes');
            all_axes = all_axes(arrayfun(@(ax) strcmp(get(ax, 'Tag'), ''), all_axes));

            
            all_pos = arrayfun(@(ax) get(ax, 'Position'), all_axes, 'UniformOutput', false);
            all_pos = vertcat(all_pos{:});

            
            y_centers = all_pos(:,2) + all_pos(:,4)/2;
            [~,~,row_ids] = uniquetol(y_centers, 1e-3);
            n_rows = max(row_ids);

            
            max_width = 0;
            ref_row_pos = [];
            ref_row_axes = [];

            for row = 1:n_rows
                idx = find(row_ids == row);
                row_axes = all_axes(idx);
                row_pos  = all_pos(idx, :);
                [~, sort_idx] = sort(row_pos(:,1));
                row_pos = row_pos(sort_idx, :);
                row_axes = row_axes(sort_idx);

                x1 = row_pos(1,1);
                x2 = row_pos(end,1) + row_pos(end,3);
                width = x2 - x1;

                if width > max_width
                    max_width = width;
                    ref_row_pos = row_pos;
                    ref_row_axes = row_axes;
                end
            end

            
            n_ref = size(ref_row_pos,1);
            subplot_width = ref_row_pos(1,3);
            gap_list = zeros(1, n_ref - 1);
            for i = 1:n_ref - 1
                left = ref_row_pos(i,1) + ref_row_pos(i,3);
                right = ref_row_pos(i+1,1);
                gap_list(i) = right - left;
            end
            gap_sum = sum(gap_list);
            if gap_sum < 1e-10
                gap_ratios = ones(1, length(gap_list)) / length(gap_list);
            else
                gap_ratios = gap_list / gap_sum;
            end

            
            usable_width = 1 - left_margin - right_margin;
            total_subplot_width = subplot_width * n_ref;
            expanded_gap_total = usable_width - total_subplot_width;
            expanded_gaps = gap_ratios * expanded_gap_total;

            
            for row = 1:n_rows
                idx = find(row_ids == row);
                row_axes = all_axes(idx);
                row_pos  = all_pos(idx, :);

                [~, sort_idx] = sort(row_pos(:,1));
                row_axes = row_axes(sort_idx);
                row_pos  = row_pos(sort_idx, :);

                n = numel(row_axes);
                x = left_margin;

                for i = 1:n
                    new_pos = row_pos(i,:);
                    new_pos(1) = x;
                    set(row_axes(i), 'Position', new_pos);
                    x = x + subplot_width;
                    if i < n
                        
                        if length(expanded_gaps) >= i
                            x = x + expanded_gaps(i);
                        else
                            x = x + expanded_gaps(end); % 
                        end
                    end
                end
            end
        end



end
    %%
    methods (Access=private)
 
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

            root = fileparts(fileparts(mfilename("fullpath")));
            addpath([root '/toolboxs/brewermap'])
            addpath([root '/toolboxs/export_fig-3.34'])
            addpath([root '/toolboxs/plotSpread'])
            addpath([root '/toolboxs/scottclowe-superbar-ce333a9/superbar'])
 
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
            % check if the 'target' in inputVars match the tarValue
            % if not, overwirte tarValue by inputVars
            % if target not exist in inputVars add it to pararms 
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
  
    end

end
