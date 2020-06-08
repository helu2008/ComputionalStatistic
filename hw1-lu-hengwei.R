# Hengwei Lu
# A99010013
# Question 1.A
load("earthquakes-2014.rda")
dat_mag <- subset(dat, magnitude >= 2)
# produce a table of counts for the number of events in each month
table(dat_mag$month)
# draw a relevant plot of month and events
hist(dat_mag$month, main = "Earthquakes in 2014", xlab = "month")
# There were more earthquakes in winter or summer. Probably earthquakes are related to temperature.

# Question 1.B
# use chi-square test, because
# Null hypothesis: earthquakes are independently distributed in every month.
chisq.test(dat_mag$month)
# p-value < 2.2e-16, so we will reject the null hypothesis. Earthquakes do not distributed evenly in every month.

# Question 2
# Yes. We can use chi-square test.
# Null hypothesis: the applicants' gender didn't have an influence on admission decision.
# male students applied
male_apply <- 8442
# male students accepted
male_accept <- 3738
# female students applied
female_apply <- 4321
# female students accepted
female_accept <- 1494
# create a contingency table
accept_table <- matrix(c(male_accept, male_apply- male_accept, 
                         female_accept, female_apply- female_accept), ncol = 2)
rownames(accept_table) <- c("accepted", "not accepted")
colnames(accept_table) <- c("male", "female")
accept_table <- as.table(accept_table)
# chi`square test
chisq.test(accept_table)
# p-value < 2.2e-16, so we will reject the null hypothesis. The applicants' gender do have an influence on admission decision.

# Question 3.A
# sum all 6 departments
sum_accept_table <- rowSums(UCBAdmissions, dims=2)
# chi-square test
# Null hypothesis: the applicants' gender didn't have an influence on admission decision.
chisq.test(sum_accept_table)
# p-value < 2.2e-16, so we will reject null hypothesis. 
# Question 3.B
plot(as.table(t(UCBAdmissions[,,1])), main = "Department A")
# In department A, more female were admitted than male..
plot(as.table(t(UCBAdmissions[,,2])), main = "Department B")
# In department B, more male were admitted than female. But they were almost the same.
plot(as.table(t(UCBAdmissions[,,3])), main = "Department C")
# In department C, more female were admitted than male. But they were almost the same.
plot(as.table(t(UCBAdmissions[,,4])), main = "Department D")
# In department D, more male were admitted than female. But they were almost the same.
plot(as.table(t(UCBAdmissions[,,5])), main = "Department E")
# In department E, more female were admitted than male. But they were almost the same.
plot(as.table(t(UCBAdmissions[,,6])), main = "Department F")
# In department F, more male were admitted than female. But they were almost the same.
# To conclude, it was not congruent because in some departments, more female were admitted compared to applied.

# Question 3.C
# Null hypothesis: Acceptance rates in all departments are independent of applicants' genders.
chisq.test(matrix(aperm(UCBAdmissions, c(2,1,3)), nrow = 2))
# p-value < 2.2e-16, so we will reject null hypothesis.