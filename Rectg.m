function r =  Rectg(xmin,ymin,xmax,ymax)
%生成矩形数据结构

r = [2 xmin xmax ymin ymin 1 0;
2 xmax xmax ymin ymax 1 0;
2 xmax xmin ymax ymax 1 0;
2 xmin xmin ymax ymin 1 0]';