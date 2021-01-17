clear all
close all
% Matlab-2 Assignment Problem 1.a

random_bits = round(rand(1,10^5));
T = 10^(-3);
s1 = zeros(1,40); s2 = zeros(1,40);

s1(1:20) = -1; s1(21:40) = 1; % bit '0'
s2(1:20) = 2; s2(21:40) = -2; % bit '1'

i = 1;
s = zeros(1, 40*length(random_bits)); % All bits are added because we need more than 5 symbols in further parts
while i <= length(random_bits)
    if random_bits(i) == 0
        s(40*(i-1) + 1: 40*i) = s1;
    else
        s(40*(i-1) + 1: 40*i) = s2;
    end
    i = i + 1;
end

t = 0: T/40: 199*T/40;
figure
plot(t, s(1:200))
title("" + random_bits(1) + random_bits(2) + random_bits(3) + random_bits(4) + random_bits(5))

% Problem 1.b

Eb = 5/2*10^(-3); %energy theoretical

low_var = 10^(0); % SNR = 0.0025 SNR_dB = -26.02 dB
med_var = 10^(-2); % SNR = 0.25 SNR_dB = -6.02  dB
high_var = 10^(-4); % SNR = 25 SNR = 13.98 dB

s_high_SNR = s + sqrt(high_var) * randn(size(s));
s_med_SNR = s + sqrt(med_var) * randn(size(s));
s_low_SNR = s + sqrt(low_var) * randn(size(s));

figure
plot(0: T/40: 199*T/40, s_high_SNR(1:200))
title('HIGH SNR')
figure
plot(0: T/40:199*T/40, s_med_SNR(1:200))
title('MEDIUM SNR')
figure
plot(0: T/40: 199*T/40, s_low_SNR(1:200))
title('LOW SNR')

% Problem 1.d

phi1 = s1 / sqrt(s1 * s1.');

s1_vector = s1 * phi1.';
s2_vector = s2 * phi1.';

correct = zeros(1, length(random_bits));
decided = zeros(1, length(random_bits));

for i = 1:10^5
    r_vector = s_med_SNR((40*(i-1) + 1: 40*i)) * phi1.';
    if norm(r_vector-s1_vector) <= norm(r_vector - s2_vector)
        decided(i) = 0;
    else
        decided(i) = 1;
    end
    if random_bits(i) == decided(i)
        correct(i) = 1;
    end
end

Pe_experimental = 1 - sum(correct) / length(correct);

disp("Probability of error: " + Pe_experimental)

% Problem 1.e

var = linspace(40000, 9.5, 10000);
Pe_exp_total = zeros(1, 10000);

for i = 1:length(var)
    total_correct = zeros(size(random_bits));
    total_decided = zeros(size(random_bits));
    received = s + sqrt(var(i)) * randn(size(s));
    for j = 1:length(random_bits)
        r_vector = received((40*(j-1) + 1: 40*j))*phi1.';
        if norm(r_vector - s1_vector) <= norm(r_vector - s2_vector)
            total_decided(j) = 0;
        else
            total_decided(j) = 1;
        end
        if random_bits(j) == total_decided(j)
            total_correct(j) = 1;
        end
    end
    Pe_exp_total(i) = 1 - sum(total_correct) / length(total_correct);
end

avg_energy = 0.5 * (s1_vector^2 + s2_vector^2);
SNR = avg_energy ./ var;
var_theo = 1/400 * var ./ avg_energy;
figure
hold on
plot(SNR, Pe_exp_total)
plot(SNR, qfunc(3 ./ (20 * sqrt(10) * sqrt(var_theo))))
title('Experimental vs. Theoretical Error')
legend('Experiment','Theory')

% % Problem 2.a
% 
% random_bits_2 = round(rand(1, 2*10^5));
% Ts = 0.2;
% d = 3;
% fc = 100;
% t = linspace(0, Ts, 180);
% si_t = zeros(16, 180);
% si = [-3*d 3*d; -3*d d; -3*d -d; -3*d -3*d; -d 3*d; -d d; -d -d; -d -3*d; d 3*d; d d; d -d; d -3*d; 3*d 3*d; 3*d d; 3*d -d; 3*d -3*d];
% 
% for i = 1:16
%         si_t(i,:) = sqrt(2/Ts)*(si(i,:)*[cos(2*pi*fc*t); -sin(2*pi*fc*t)]);
% end
% 
% random_4_bits = strings(1,length(random_bits_2)/4);
% resulted_bits = ["0000" "0001"  "0010" "0011" "0100" "0101" "0110" "0111" "1000" "1001" "1010" "1011" "1100" "1101" "1110" "1111"];
% resulted_gray_bits = ["0000" "0001" "0011" "0010" "0100" "0101" "0111" "0110" "1100" "1101" "1111" "1110" "1000" "1001" "1011" "1010"];
% for i = 1:length(random_bits_2)/4
%     random_4_bits(i) = "" + random_bits_2(4*i-3) + random_bits_2(4*i-2) + random_bits_2(4*i-1) + random_bits_2(4*i);
% end
% 
% for i = 1:length(random_bits_2)/4
%     for j = 1:16
%     if random_4_bits(i) == resulted_bits(j)
%         s((180*(i-1) + 1): 180*i) = si_t(j,:);
%     end
%     end
% end
% 
% for i = 1:length(random_bits_2)/4
%     for j = 1:16
%     if random_4_bits(i) == resulted_gray_bits(j)
%         s_gray((180*(i-1) + 1): 180*i) = si_t(j,:);
%     end
%     end
% end
% figure
% plot( 0: Ts/180: Ts*(5-1/180), s(1:900) )
% title("" + random_4_bits(1) + " " + random_4_bits(2) + " " + random_4_bits(3) + " " + random_4_bits(4) + " " + random_4_bits(5))
% 
% % Problem 2.b
% 
% high_var = 5;
% 
% s_high_SNR = s + sqrt(high_var)*randn(size(s));
% s_high_SNR_gray = s_gray + sqrt(high_var)*randn(size(s));
% 
% figure
% plot( 0: Ts/180: Ts*(5-1/180),s_high_SNR(1:900) )
% title('Noisy Signal')
% 
% phi1_2 = cos(2*pi*fc*t)/(sqrt(cos(2*pi*fc*t)*cos(2*pi*fc*t).'));
% phi2_2 = sin(2*pi*fc*t)/(sqrt(sin(2*pi*fc*t)*sin(2*pi*fc*t).'));
% 
% constellation = si_t * phi1_2.';
% constellation_2 = si_t * phi2_2.';
% 
% r1 = zeros(1,length(s)/180);
% r2 = zeros(1,length(s)/180);
% for i = 1:length(s)/180
%     r1(i) = s_high_SNR(180*(i-1) + 1: 180*i) * phi1_2.';
%     r2(i) = s_high_SNR(180*(i-1) + 1: 180*i) * phi2_2.';
% end
% 
% r1_gray = zeros(1, length(s)/180);
% r2_gray = zeros(1, length(s)/180);
% for i = 1:length(s)/180
%     r1_gray(i) = s_high_SNR_gray(180*(i-1) + 1: 180*i) * phi1_2.';
%     r2_gray(i) = s_high_SNR_gray(180*(i-1) + 1: 180*i) * phi2_2.';
% end
% 
% figure
% hold on
% for i = 1:100
%     plot(r1(i), r2(i), '*r')
% end
% for i = 1:16
%     plot(constellation(i), constellation_2(i), '*b')
% end
% title('r1 r2 with Constellation Points')
% 
% % Problem 2.c
% 
% r = [r1;r2].';
% s_vector = [constellation constellation_2];
% norm_arr = zeros(1, 16);
% estimated_message_4b = strings(1, length(random_4_bits));
% estimated_bits = zeros(1, length(random_bits_2));
% total_correct_4 = zeros(1, length(random_bits_2)/4);
% total_correct_bits = zeros(1, length(random_bits_2));
% 
% for i = 1:length(random_bits_2)/4
%     for j = 1:16
%         norm_arr(j) = norm( r(i,:) - s_vector(j,:) );
%     end
%     estimate_index = find(norm_arr == min(norm_arr));
%     estimated_message_4b(i) = resulted_bits(estimate_index);
%     if random_4_bits(i) == estimated_message_4b(i)
%             total_correct_4(i) = 1;
%     end
% end
% 
% chr = char(estimated_message_4b);
% for i = 1:length(random_bits_2)
%     temp1 = str2double(chr(i));
%     estimated_bits(i) = temp1;
%     if random_bits_2(i) == estimated_bits(i)
%             total_correct_bits(i) = 1;
%     end
% end
% 
% % gray coding
% r_gray = [r1_gray;r2_gray].';
% norm_arr_gray = zeros(1, 16);
% estimated_message_4b_gray = strings(1, length(random_4_bits));
% estimated_bits_gray = zeros(1, length(random_bits_2));
% total_correct_bits_gray = zeros(1, length(random_bits_2));
% 
% for i = 1:length(random_bits_2)/4
%     for j = 1:16
%         norm_arr_gray(j) = norm( r_gray(i,:) - s_vector(j,:) );
%     end
%     estimate_index_gray = find(norm_arr_gray == min(norm_arr_gray));
%     estimated_message_4b_gray(i) = resulted_gray_bits(estimate_index_gray);
% end
% 
% chr_gray = char(estimated_message_4b_gray);
% for i = 1:length(random_bits_2)
%     temp1 = str2double(chr_gray(i));
%     estimated_bits_gray(i) = temp1;
%     if random_bits_2(i) == estimated_bits_gray(i)
%             total_correct_bits_gray(i) = 1;
%     end
% end
% 
% Pe_exp_bit = 1 - sum(total_correct_bits) / length(random_bits_2);
% Pe_exp_symbol = 1 - sum(total_correct_4) / length(random_4_bits);
% disp("Statistical probability of symbol error is " + Pe_exp_symbol)
% disp("Statistical probability of bit error: " + Pe_exp_bit)
% Pe_exp_gray_bit = 1 - sum(total_correct_bits_gray) / length(random_bits_2);
% disp("Statistical probability of gray bit error: " + Pe_exp_gray_bit)
% 
% % Problem 2.d
% 
% var = linspace(12100, 670, 10);
% Pe_exp_symbol_total = zeros(1, 10);
% estimated_m = strings(1, length(random_4_bits));
% 
% for j = 1:length(var)
%     estimate_check = zeros(1,length(random_bits_2)/4);
%     received = s + sqrt(var(j))*randn(size(s));
%     for i = 1:length(random_bits_2)/4
%         received_r = [received((180*(i-1) + 1: 180*i)) * phi1_2.' received((180*(i-1) + 1: 180*i)) * phi2_2.'];
%         norm_arr = zeros(1,16);
%         for k = 1:16
%             norm_arr(k) = norm(received_r - s_vector(k,:));
%         end
%         estimate_ind = find(norm_arr == min(norm_arr));
%     estimated_m(i) = resulted_bits(estimate_ind);
%     if random_4_bits(i) == estimated_m(i)
%             estimate_check(i) = 1;
%     end
%     end
%     Pe_exp_symbol_total(j) = 1 - sum(estimate_check)/length(estimate_check);
% end
% 
% avg_energy = sum(constellation.^2 + constellation_2.^2) / 16;
% avg_energy_theo = 90;
% SNR = avg_energy./var;
% sigma_theoretical = sqrt(avg_energy_theo * var / avg_energy);
% figure
% hold on
% plot(SNR,Pe_exp_symbol_total)
% plot(SNR,3 * qfunc(3 ./ (sigma_theoretical)) - 9/4 * (qfunc(3 ./ (sigma_theoretical))).^2);
% title('Experimental Error vs. Theoretical Error')
% legend('Experimental','Theoretical')