function [state,options,optchanged] = gaoutfun(options,state,flag)

optchanged = false;

persistent GA_POPhistory GA_scorehistory GA_ERRhistory

switch flag
    case 'init'
        clear GA_POPhistory
        clear GA_scorehistory
        clear GA_ERRhistory
        
    case 'iter'
        % GA_POPhistory(:,:,state.Generation) = state.Population;
        GA_scorehistory(:,state.Generation) = state.Score;
        GA_ERRhistory(1:state.Generation) = state.Best;
        
    case 'done'
        % save('tempResultsFile.mat','GA_POPhistory','GA_scorehistory','GA_ERRhistory')
        save('tempResultsFile.mat','GA_scorehistory','GA_ERRhistory')
        
end