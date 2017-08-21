function [res] = u( a, x, y )
    res = -a * y / ( x^2 + y^2 );
end