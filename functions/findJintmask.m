function [Jint] = findJintmask(~, mesh, el, Jint)

nodes = reshape(1:length(mesh.UDIC),mesh.winDIC(2),mesh.winDIC(1));

Jintnodes = ones(size(nodes));

for i=1:size(Jint.mask,3)
    Jintnodes(min(Jint.mask(1,:,i)):max(Jint.mask(1,:,i)),...
        min(Jint.mask(2,:,i)):max(Jint.mask(2,:,i)))=0;
end

Jintnodes = Jintnodes.*nodes;
Jintnodes = unique(Jintnodes);
Jintnodes(Jintnodes==0)=[];

elnum =1;
for c = 1:length(el.n)
    if(sum(ismember(el.n(c,:),Jintnodes))==4)
        JDelem(elnum) = c;
        elnum = elnum+1;
    end
end

Jint.unmasked = JDelem';
