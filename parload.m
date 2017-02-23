% load vartiables 'normedscore10', 'x', and 'y' inside a parfor loop
function [data, x, y] = parload(filename)
    data1 = load(filename, 'normedscore10');
    data = data1.normedscore10;
    x1 = load(filename, 'x');
    x = x1.x;
    y1 = load(filename, 'y');
    y = y1.y;
end