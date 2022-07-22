% Converts edge list to adjacency matrix
% INPUTS: edgelist: mx3
% OUTPUTS: adjacency matrix nxn
% Note: information about nodes is lost: indices only (i1,...in) remain
% GB, Last updated: October 6, 2009

% function adj=edgeL2adj(el)
% 
% nodes=sort(unique([el(:,1) el(:,2)])); % get all nodes, sorted
% adj=zeros(numel(nodes));   % initialize adjacency matrix
% 
% % across all edges
% for i=1:size(el,1); adj(find(nodes==el(i,1)),find(nodes==el(i,2)))=el(i,3);
% adj(find(nodes==el(i,2)),find(nodes==el(i,1)))=el(i,3);
% end


function adj=edgeL2adj0(el)

nodes=sort(unique([el(:,1) el(:,2)])); % get all nodes, sorted
adj=zeros(numel(nodes));   % initialize adjacency matrix

% across all edges
for i=1:size(el,1); 
    adj(find(nodes+1==el(i,1)+1),find(nodes+1==el(i,2)+1))=1;
    adj(find(nodes+1==el(i,2)+1),find(nodes+1==el(i,1)+1))=1;
end