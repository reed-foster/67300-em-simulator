function f = color_spy(M)
  % plots SPY of matrix but with colorcoded entries for
  a = M(:,:);
  a(a > 0) = 1;
  a(a < 0) = -1;
  f = figure;
  colormap(brewermap([],"RdBu"));
  h = pcolor(a);
  set(h, 'EdgeColor', 'none');
  axis ij
  axis square;
end
