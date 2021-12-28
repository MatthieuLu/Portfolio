function status = func_collision(~,S,flag)

if size(S)==[0,0]
    return
end

Z=S(1,:)';

for i=1:size(Z,1)
    
    z=Z(i,1);
    
%     if strcmp(flag,'done')
%         return
%     end

    if z<0
        fprintf('Collision with the earth!!!\n\n')
        status=1;
        flag = 'done';
        return
    else
        status=0;
    end
    
end