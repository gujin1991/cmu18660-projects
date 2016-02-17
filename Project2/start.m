tic;
imgIn = imgRead('fishing_boat.bmp');
%imgIn = imgRead('lena.bmp');
imgOut = imgRecover(imgIn, 8, 10);
imgShow(medfilt2(imgOut, [3 3]));
error = mean(mean((imgOut - imgIn) .^2)); %initial errr
medError = mean(mean((medfilt2(imgOut, [3 3]) - imgIn) .^2)); %error after medium filter
toc;