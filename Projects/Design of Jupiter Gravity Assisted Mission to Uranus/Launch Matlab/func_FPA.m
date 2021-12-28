function status = func_FPA(~,S,flag)

if size(S)==[0,0]
    return
end

FPA=S(3,:)';

for i=1:size(FPA,1)
    
    fpa=FPA(i,1);
    
    if fpa<0
        status=1;
        flag = 'done';
        return
    else
        status=0;
    end
    
end