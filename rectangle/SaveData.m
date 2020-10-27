function SaveData(x, y)
% intpu: x and y - arrays for points
dlmwrite('200_anis_x_y_z.txt', [x' y'],'delimiter','\t');
save 'rect_pdf_output';
end