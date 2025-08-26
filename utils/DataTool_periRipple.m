classdef DataTool_periRipple
    %DataTool_periRipple Summary of this class goes here
    %   Detailed explanation goes here

    properties
        D_signal
        D_chanTbl
        D_ROI
        D_chanid
        D_times 
        D_subject
        D_general
        D_general_sublevel
        D_gWholeTask
    end

    methods
        function obj = DataTool_periRipple(inputArg1,inputArg2)
            %DataTool_periRipple Construct an instance of this class
            %   Detailed explanation goes here
            obj.D_signal   = evalin('base','ROI_DATA');
            obj.D_chanTbl  = evalin('base','chanTbl');
            obj.D_chanid = evalin('base','ROI_CHANNELID');
            obj.D_ROI      = evalin('base', 'ROI');
            obj.D_times      = evalin('base', 'times');
            obj.D_subject  =  evalin('base', 'subjects');
            obj.D_subject  = strcat('sub-',regexprep(obj.D_subject.list,'.*sub-',''));
            obj.D_general  = evalin('base','EEG_ALL_generalAvg');
            obj.D_gWholeTask  = evalin('base','EEG_ALL_generalAvg');
            obj.D_general_sublevel = f_extract_by_ROI(obj,obj.D_general);
        end
        function output = f_clean_EC_from_neocortex(obj)
            warning('Entorhinal cortex channel will be exclude from Yeo-7')
            chanid_EC = obj.D_chanid.Entorhinal{1};
            rois = setdiff(obj.D_signal.Properties.VariableNames,'Entorhinal');
            for i = 1:numel(rois)
                iroi = rois{i};
                ind_select = ~ismember(obj.D_chanid.(iroi){1},chanid_EC);
                obj.D_signal.(iroi){1} = obj.D_signal.(iroi){1}(ind_select,:);
                obj.D_chanid.(iroi){1} = obj.D_chanid.(iroi){1}(ind_select);
                obj.D_chanTbl.(iroi) = obj.D_chanTbl.(iroi)&(~obj.D_chanTbl.Entorhinal);
                
                ind_logic = obj.D_chanTbl.(iroi);
                [m,gname] = grpstats(obj.D_general(ind_logic), ...
                    obj.D_chanTbl.subject(ind_logic),{@nanmean,'gname'});
                obj.D_general_sublevel.(iroi) = nan(size(obj.D_subject) );
                [~,ia,ib] = intersect(obj.D_subject,gname,'stable');
                obj.D_general_sublevel.(iroi)(ia) = m(ib);
            end
            output = obj;
            

            
        end
        function outputArg = f_subLevel_timeCourse(obj)
                roinames = obj.D_signal.Properties.VariableNames ;
                amp ={}; outputArg=table;
                for iroi = 1:numel(roinames)
                    roi_name = roinames{iroi};
                    d = obj.D_signal.(roi_name){1};
         
                    roi_ind    = obj.D_chanTbl.(roi_name);              
                    channel   = obj.D_chanTbl.channelid(roi_ind) ;

                    subject = regexprep(channel,'.*_','');
                    [gm,gname] = groupsummary(d,subject,@mean);
                    [~,ia,ib] = intersect(obj.D_subject,gname,'stable');
                    tmp = nan(numel(obj.D_subject),size(d,2));
                    tmp(ia,:)= gm(ib,:);
                    amp{iroi} = tmp;
                    outputArg.(roi_name)=tmp;
                end
%             outputArg = cell2table(amp,'VariableNames',roinames);
              
        end


        function outputArg = f_extractByWindow_and_ROI(obj,windows)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if ~exist('windows','var')||isempty(windows)
                windows = [ -0.15 0.15; -0.25 0.25; -0.5 0.5;];
            end

            outputArg = table; 
            roinames = obj.D_signal.Properties.VariableNames ;
            for iwind = 1: size(windows,1)
                DW = struct;
                for iroi = 1:numel(roinames)
                    roi_name = roinames{iroi};
                    d = obj.D_signal.(roi_name){1};
         
                    roi_ind    = obj.D_chanTbl.(roi_name);              
                    channel   = obj.D_chanTbl.channelid(roi_ind) ;

                    bin = windows(iwind,:)*1000;
                    ind = obj.D_times>bin(1)&obj.D_times<bin(2);
   
                    amp = mean(d(:,ind),2);
 
                    amp = table(amp,channel);

                    DW.(roi_name)= amp;
  
                end

                outputArg.(sprintf('wind_%.02f~%.02f',windows(iwind,:)))=DW;
            end
        end

        function varargout = subjectLevel(obj,sortROIbyAvg,varargin)
            %
            obj = obj.f_clean_EC_from_neocortex();
            % defaults
            if ~exist('sortROIbyAvg','var'),sortROIbyAvg =true;end
            data = obj.f_extractByWindow_and_ROI(varargin{:});
            windows = data.Properties.VariableNames;
            rois = fieldnames(data{1,1});
            outputArg = table;outputArg1 = table;outputArg2 = table;
            for iwind = 1: numel(windows)
                d = data.(windows{iwind});
                amp  = nan(numel(obj.D_subject),numel(rois));
                for iroi = 1:numel(rois)
                    dd = d.(rois{iroi});
                    subject = regexprep(dd.channel,'.*_','');
                    [gm,gname] = grpstats(dd.amp,subject,{'mean','gname'});
                    [~,ia,ib] = intersect(obj.D_subject,gname,'stable');
                    amp(ia,iroi) = gm(ib);
                end            
                amp=array2table(amp,'VariableNames',rois);
                if sortROIbyAvg
                    amp_mean = nanmean(amp.Variables,1);
                    [~,order] = sort(amp_mean,'descend');
                    amp = amp(:,order);
                end
%                 amp.subject = obj.D_subject;
                outputArg.(windows{iwind})=amp;
                if nargout>1
                    ind = contains(amp.Properties.VariableNames,'Y7_');
                    outputArg1.(windows{iwind})=amp(:,ind);
                    outputArg2.(windows{iwind})=amp(:,~ind);
                end
            end
            outputArg.subject = obj.D_subject;
            if nargout>1
                outputArg1.subject = obj.D_subject;
                outputArg2.subject = obj.D_subject;
                varargout = {outputArg1,outputArg2};
            else
                varargout  = {outputArg};
            end
            varargout{end+1}=order;

        end
        function varargout = channelLevel(obj,sortROIbyAvg,varargin)
            % defaults
            if ~exist('sortROIbyAvg','var'),sortROIbyAvg =true;end

            data = obj.f_extractByWindow_and_ROI(varargin{:});
            windows = data.Properties.VariableNames;
            rois = fieldnames(data{1,1});
            outputArg = table;outputArg1 = table;outputArg2 = table;
            for iwind = 1: numel(windows)
                d = data.(windows{iwind});
                channel_all =(obj.D_chanTbl.channelid);
                amp  = nan(numel(channel_all),numel(rois));
                for iroi = 1:numel(rois)
                    dd = d.(rois{iroi});
 
                    [~,ia,ib] = intersect(channel_all,dd.channel,'stable');
                    amp(ia,iroi) =dd.amp(ib);
                end            
                amp=array2table(amp,'VariableNames',rois);
                if sortROIbyAvg
                    amp_mean = nanmean(amp.Variables,1);
                    [~,order] = sort(amp_mean,'descend');
                    amp = amp(:,order);
                end
%                 amp.subject = obj.D_subject;
                outputArg.(windows{iwind})=amp;
                if nargout>1
                    ind = contains(amp.Properties.VariableNames,'Y7_');
                    outputArg1.(windows{iwind})=amp(:,ind);
                    outputArg2.(windows{iwind})=amp(:,~ind);
                end
            end
            outputArg.channel = channel_all;
            if nargout>1
                outputArg1.channel =channel_all;
                outputArg2.channel =channel_all;
                varargout = {outputArg1,outputArg2};
            else
                varargout  = {outputArg};
            end
            

        end        
        function myCM = f_colormap(obj)
            myCM.vmPFC=[.5 .5 .5];
        	myCM.dmPFC=[.5 .5 .5];
            myCM.mPFC =[.5 .5 .5];
        	myCM.vmPFC_Large=[.5 .5 .5];
            myCM.dmPFC_Large=[.5 .5 .5];
            myCM.mPFC_Large =[.5 .5 .5];
            myCM.DMNPCC  =[.5 .5 .5];            
            myCM.DMNmPFC =[.5 .5 .5];
            myCM.DMNmPFC =[.5 .5 .5];
        	myCM.Entorhinal =[.5 .5 .5];
            myCM.Y7_Visual  =[0.46875    0.070312     0.52344];
            myCM.Y7_Default =[0.80078     0.24219     0.30469];
            myCM.Y7_DorsalAttention  =[      0     0.46094    0.054688];
            myCM.Y7_Frontoparietal  =[0.89844     0.57812     0.13281];
            myCM.Y7_VentralAttention  =[0.76562     0.22656     0.97656];
            myCM.Y7_Limbic  =[0.85938     0.96875     0.64062];
            myCM.Y7_Somatomotor  =[0.27344     0.50781     0.70312];
  % 'DMNmPFC' 'DMNPCC' 'DMNDLPFC' 'DMNLOFC' 'DMNTPJ' 'DMNTL'          
            myCM = struct2table(myCM);
        end        
        function D_general_sublevel = f_extract_by_ROI(obj,d)
            D_general_sublevel = table;
            for iroi = 1:numel(obj.D_ROI.names)
                iroiname = obj.D_ROI.names{iroi};
                
                chanind = obj.D_chanTbl.(iroiname);
                subject = obj.D_chanTbl.subject(chanind);
                

                [dd,gname] =  grpstats(d(chanind,:),subject,{'mean','gname'});
                [~,ia,ib] = intersect(obj.D_subject,gname,'stable');
                amp = nan(size(obj.D_subject,1),size(d,2));
                amp(ia,:) = dd(ib,:);
                D_general_sublevel.(iroiname) = amp;        
            end            
        end        
    end

    methods (Access=private)
        function obj = f_add_mPFC_inDMN(obj)
            

        end
        
    end
end