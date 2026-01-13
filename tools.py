import matplotlib.pyplot as plt

def draw(data_series, xlabel='', ylabel='', title='', xlim=None, ylim=None):
    fig, ax = plt.subplots(figsize=(10, 6))
    for series in data_series:
        x = series.get("x", [])
        y = series.get("y", [])
        style = series.get("style", 'b-')
        label = series.get("label", '')
        ax.plot(x, y, style, linewidth=2, label=label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.savefig('out/pressure_velocity_curve.png', dpi=300)

    plt.tight_layout()
    plt.show()
