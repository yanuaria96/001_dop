function SaveData(x, y)
dlmwrite('200_anis_x_y_z.txt', [x' y'],'delimiter','\t');
save 'rect_pdf_output';
end