function [ ip ] = permutation( i,j,k,pl,np )
%PERMUTATION Summary of this function goes here
%   Stolen from ATHENA to calculate point weights. No clue how it works
    ip=-1;
    
  for l=1:np
    for m=1:3
      if (i == pl(l,m))
        for n=1:3
          if(n ~= m)
            if (j == pl(l,n))
              for o=1:3
                if((o ~= m) && (o ~= n))
                  if (k == pl(l,o))
                    ip = l;
                  end
                end
              end
            end
          end
        end
      end
    end
  end

end

