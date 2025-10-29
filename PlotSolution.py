#!/usr/bin/env python
# coding: utf-8

# In[6]:


import matplotlib.pyplot as plt
import numpy as np
def plot(customers, facilities):
    assignments = {}
    for facility in facilities:
        assignments[facility] = []

    for customer in customers:
        distances = [np.hypot(customer.latitude - facility.latitude,
                              customer.longitude - facility.longitude) for facility in facilities]
        nearest = facilities[np.argmin(distances)]
        assignments[nearest].append(customer)
    
    plt.figure(figsize=(6, 6))
    plt.title("Customer-Facility Connections")
    
    # Plot facilities
    for facility in facilities:
        plt.scatter(facility.latitude, facility.longitude, color='red', s=100, marker='^')
        plt.text(facility.latitude + 0.1, facility.longitude, facility.name, fontsize=10, color='black')
    
    # Plot customers and connections
    for facility, customers in assignments.items():
        for customer in customers:
            plt.scatter(customer.latitude, customer.longitude, color='blue', s=30, linewidths=0.1)
            plt.plot([facility.latitude, customer.latitude], 
                     [facility.longitude, customer.longitude], color='blue', alpha=0.3)
    
    plt.xlabel("Latitude")
    plt.ylabel("Longitude")
    plt.grid(True, alpha=0.3)
    plt.show()

