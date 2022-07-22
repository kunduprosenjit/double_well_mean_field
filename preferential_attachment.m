
function el = pref_attach(n,m)

ver = 2;

el=[1 2 1; 2 1 1];  % start with an edge


while ver < n
    ver=ver+1;  % add new vertex

    if m>=ver
      for node=1:ver-1
        el = [el; node ver 1];
        el = [el; ver node 1];
      end
      continue
    end
    
    deg=[];               % compute nodal degrees for this iteration
    for v=1:ver 
        deg=[deg; v numel(find(el(:,1)==v))];
    end
    deg=sortrows(deg);
    
    % add m edges
    r=randsample(deg(:,1),m,'true',deg(:,2)/max(deg(:,2)));
    while not(length(unique(r))==length(r))
      r=randsample(deg(:,1),m,'true',deg(:,2)/max(deg(:,2)));
    end

    for node=1:length(r)
      el=[el; r(node) ver 1];
      el=[el; ver r(node) 1];      
    end      
    
end