% edge list to adjacency matrix
% INPUTS: edge list: nx3
% OUTPUTS: adjacency matrix nxn

%egl=edge list

function adj=edgeL2adj(egl)

nodes=sort(unique([egl(:,1) egl(:,2)])); % get all nodes, sorted
adj=zeros(numel(nodes));   % initialize adjacency matrix

% if node indices starts from zero.
for i=1:size(egl,1); 
    adj(find(nodes+1==egl(i,1)+1),find(nodes+1==egl(i,2)+1))=1; 
    adj(find(nodes+1==egl(i,2)+1),find(nodes+1==egl(i,1)+1))=1;
end

% % if node indices starts from 1.
% for i=1:size(egl,1); 
%     adj(find(nodes==egl(i,1)),find(nodes+1==egl(i,2)))=1; 
%     adj(find(nodes==egl(i,2)),find(nodes+1==egl(i,1)))=1;
% end
