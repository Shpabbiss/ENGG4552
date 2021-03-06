import numpy as np
import matplotlib.pyplot as plt
scale = 1
diameter = 300
radius = int((diameter*scale)/2)
n = 90
width = 100*scale
Nx = Ny = int(radius*2+np.sqrt(radius**2 + width/2*np.sin(((180/n)-90)))+1)
pixel_array = np.zeros((Nx, Ny))
# pixel_array[501, 501] = 10
print(np.shape(pixel_array))
theta = 180 / n
centre_point = (int(Nx/2), int(Ny/2))

def calc_side(n):
    return 2 * radius * np.sin(np.pi / n)


# def cacl_rect_strip(theta):
#     centre_point = (501,501)
#     circle = pixel_array.copy()
#     for r in range(radius):
#         x_loc = int(centre_point[0]+r*np.sin(theta))
#         y_loc = int(centre_point[1]+r*np.cos(theta))
#         circle[x_loc, y_loc] +=1
#     return circle

def cacl_rect_strip(theta):

    circle = pixel_array.copy()

    if width % 2 != 0:
        print("Not a nice width")

    for w in range(width):
        print
        for r in range(radius):
            x_loc = int(centre_point[0] + r * np.sin(theta) - w / 2 * np.sin(theta + (np.pi / 2)))
            y_loc = int(centre_point[1] + r * np.cos(theta) - w / 2 * np.cos(theta + (np.pi / 2)))
            circle[x_loc, y_loc] = 1
            x_loc = int(centre_point[0] + r * np.sin(theta) - w / 2 * np.sin(theta - (np.pi / 2)))
            y_loc = int(centre_point[1] + r * np.cos(theta) - w / 2 * np.cos(theta - (np.pi / 2)))

            circle[x_loc, y_loc] = 1
            # x_loc = int(centre_point[0] + r * np.sin(theta + np.pi) - w / 2 * np.sin(theta + (np.pi / 2)))
            # y_loc = int(centre_point[1] + r * np.cos(theta + np.pi) - w / 2 * np.cos(theta + (np.pi / 2)))
            # circle[x_loc, y_loc] = 1
            # x_loc = int(centre_point[0] + r * np.sin(theta + np.pi) - w / 2 * np.sin(theta - (np.pi / 2)))
            # y_loc = int(centre_point[1] + r * np.cos(theta + np.pi) - w / 2 * np.cos(theta - (np.pi / 2)))
            # circle[x_loc, y_loc] = 1
    return circle


circ = pixel_array

for t in range(0,360, int(theta)):
    print(t)
    t = np.radians(t)
    temp = cacl_rect_strip(t)
    circ = np.add(circ,temp)
    # plt.imshow(cacl_rect_strip(t), cmap='hot', interpolation='none')
    # plt.colorbar()
    # plt.show()

plt.imshow(circ, origin='lower', cmap='hot', interpolation='none', vmax=n)

plt.colorbar()
nx = circ.shape[0]
no_labels = 7 # how many labels to see on axis x
step_x = int(nx / (no_labels - 1)) # step between consecutive labels
x_positions = np.arange(0,nx,step_x) # pixel count at label position
print(x_positions)
x_labels = x_positions/scale # labels you want to see
print(x_labels)
plt.xticks(x_positions, x_labels)
plt.yticks(x_positions, x_labels)
plt.show()

plt.plot(np.arange(0,Nx), centre_point[0])
plt.show()