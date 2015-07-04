function ber_x = Bernouli (x)
    if(x > 0.01)
         ber_x = x*exp(-x)/(1-exp(-x));
    elseif((x < 0) && (abs(x) > 0.01))
         ber_x = 1.0;
    else
         temp = 1.0;
         summ = temp;
         i = 0;
         while (1)
            i = i+1;
            temp = temp * x/(i+1);
            if(summ+temp == summ)
                summ = summ + temp;
                break;
            end;
            summ = summ + temp;
         end;
         ber_x = 1/summ;
    end;