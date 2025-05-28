# Calculation of Gating Probabilities for Clinical Decision Making

This interactive Shiny application helps clinical development teams evaluate early-phase efficacy data using response-based gating rules for response rates.

---

## Dynamic Input Controls

Adjust:

- **Target ORR thresholds** for Go and Stop decisions  
- **Probability thresholds** for decision confidence  
- **Sample size** and number of observed responders  
- **True response rates** for scenario evaluation  

---

## Decision Table Output

Based on the chosen design inputs and observed data, the app displays a detailed **gating decision table**, clearly marking:

- **GO**: Proceed to next phase  
- **STOP**: Terminate  
- **EVALUATE**: Consider further assessment  

---

## Real-time Probability Output

For any assumed **true ORR**, the app computes:

- Probability of **incorrect GO** decision (false positive)  
- Probability of **correct STOP** decision (true negative)  
- Probability that **further evaluation** is required (indeterminate outcome)  

---

## Interactive Visualizations

Includes two intuitive plots:

- **Cumulative Decision Probability Plot**: Shows combined probabilities for Stop, Evaluate, and Go across all response rates  
- **Decision Probability Plot**: Depicts individual curves for Go, Evaluate, and Stop probabilities with user-defined vertical reference line  
