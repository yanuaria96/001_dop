function ChangeSharp(iter)

global SHARP

if iter >400
    if rem(iter,20) == 0
        SHARP=SHARP+1;
    end
end

end