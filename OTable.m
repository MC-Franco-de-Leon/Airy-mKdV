
function [output]=OTable(m)
% M=[.1 .02 ;2 4 ];
% c={'iter 1','$$||\theta_{dt}-\theta_{dt/2}||_{l_2}$$'};
% r={'row 1','row 2'};
% makeHtmlTable(M,[],r,c)
% 
% writelatextable( [], M, '|c|c|' , [], [] );


    % Make the filename
    fn = '.\FractaFractorus';

    % Make the input data
   % m = [(10:-1:1)' rand([10,2])];
    % The table will have a column of integers, decimal values, and
    %  floats expressed in scientific notation
    mf = {'|||c%d||','c%g|','c%g|','c%g|','c%g|||'};

    % The first column is altitudes spans 2 rows
    h(1).text='$t_0$';    h(1).format='|||c||';  
    h(1).rows=[1 2];              h(1).cols=[1]; 
    % The two following columns share a heading called 'Atmospheric Data'
    h(2).text='$Log_2(\frac{||\theta_{dt}-\theta_{dt/2}||_{l_2}}{||\theta_{dt/2}-\theta_{dt/4}||_{l_2}})$'; h(2).format='c|||'; 
    h(2).rows=[1];                h(2).cols=[2]; 
    % The last column has a subheading called '$NO_2$'
    h(3).text='$Log_2(\frac{||\theta_{dt/2}-\theta_{dt/4}||_{l_2}}{||\theta_{dt/4}-\theta_{dt/4}||_{l_2}})$';           h(3).format='c|'; 
    h(3).rows=[1];                h(3).cols=[3];
    
    h(4).text='$Log_2(\frac{||\theta_{dt/4}-\theta_{dt/8}||_{l_2}}{||\theta_{dt/8}-\theta_{dt/16}||_{l_2}})$';           h(4).format='c|'; 
    h(4).rows=[1];                h(4).cols=[4];
    
    h(5).text='$Log_2(\frac{||\theta_{dt/8}-\theta_{dt/16}||_{l_2}}{||\theta_{dt/16}-\theta_{dt/21}||_{l_2}})$';           h(5).format='c|'; 
    h(5).rows=[1];                h(5).cols=[5];

    
  
    % make a caption and label
    c = 'Time Order of Convergence';
    l = 'table:atm_prof';

    output = writelatextable( fn, m, mf, h, c, l );