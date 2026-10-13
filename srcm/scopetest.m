function scopetest
x = 1;
child_assign;
fprintf('after nested assign: %d\n', x);
try
    child_arg(5);
    fprintf('arg-shadow ok, J=%d\n', J);
catch me
    fprintf('arg-shadow ERROR: %s\n', me.message);
end
    function child_assign
        x = 2;
    end
    function child_arg(J)
        J = J * 2;
    end
end
