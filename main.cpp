#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QWidget>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QFrame>
#include <QtCore/Qt>
#include <cmath>
#include <array>
#include <vector>
#include <functional>
#include <algorithm>
#include <limits>

// Constants
const double Rgas = 0.0831446261815324; // L bar /K /mol, SI 2019
const double Vguess = 0.01; // for now we're not using an initial guess -
// bracketing values are chosen in the findRoot function

// EOS coefficients
const std::array<double, 15> c = {
    2.95177298930e-2, -6.33756452413e+3, -2.75265428882e+5, 1.29128089283e-3, -1.45797416153e+2,
    7.65938947237e+4, 2.58661493537e-6, 0.52126532146e+0, -1.39839523753e+2, -2.36335007175e-8,
    5.35026383543e-3, -0.27110649951e+0, 2.50387836486e+4, 0.73226726041e+0, 1.54833359970e-2
};

// Lenard-Jones parameters
const std::array<double, 7> eps = {154.0, 510.0, 235.0, 31.2, 105.6, 124.5, 246.1};
const std::array<double, 7> sig = {3.691, 2.88, 3.79, 2.93, 3.66, 3.36, 4.35};

// Mixing parameters
std::array<std::array<double, 7>, 7> K1;
std::array<std::array<double, 7>, 7> K2;

const std::array<QString, 7> species = {
    "CH<sub>4</sub>", "H<sub>2</sub>O", "CO<sub>2</sub>", "H<sub>2</sub>",
    "CO", "O<sub>2</sub>", "C<sub>2</sub>H<sub>6</sub>"
};

// Initialize mixing parameters
void initializeMixingParameters() {
    // Initialize K1 and K2 to ones
    for (int i = 0; i < 7; ++i) {
        for (int j = 0; j < 7; ++j) {
            K1[i][j] = 1.0;
            K2[i][j] = 1.0;
        }
    }

    // Set specific values
    K1[1][2] = K1[2][1] = 0.85; // H2O-CO2
    K1[0][1] = K1[1][0] = 0.8;  // CH4-H2O
    K2[1][2] = K2[2][1] = 1.02; // H2O-CO2
}

// Brent's method for root finding
double brentMethod(std::function<double(double)> func, double a, double b,
                   double tolerance = 1e-12, int maxIterations = 1000) {
    double fa = func(a);
    double fb = func(b);

    if (fa * fb > 0) {
        // Try to find a bracketing interval by expanding
        double step = std::abs(b - a);
        for (int i = 0; i < 10; ++i) {
            a -= step;
            fa = func(a);
            if (fa * fb < 0) break;

            b += step;
            fb = func(b);
            if (fa * fb < 0) break;

            step *= 2;
        }

        if (fa * fb > 0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    if (std::abs(fa) < std::abs(fb)) {
        std::swap(a, b);
        std::swap(fa, fb);
    }

    double c = a;
    double fc = fa;
    bool mflag = true;
    double d = 0;

    for (int i = 0; i < maxIterations; ++i) {
        if (std::abs(fb) < tolerance || std::abs(b - a) < tolerance) {
            return b;
        }

        double s;
        if (fa != fc && fb != fc) {
            // Inverse quadratic interpolation
            s = a * fb * fc / ((fa - fb) * (fa - fc)) +
                b * fa * fc / ((fb - fa) * (fb - fc)) +
                c * fa * fb / ((fc - fa) * (fc - fb));
        } else {
            // Secant method
            s = b - fb * (b - a) / (fb - fa);
        }

        // Check conditions for bisection
        double tmp1 = (3 * a + b) / 4;
        double tmp2 = b;
        if (tmp1 > tmp2) std::swap(tmp1, tmp2);

        bool condition1 = !(s >= tmp1 && s <= tmp2);
        bool condition2 = mflag && std::abs(s - b) >= std::abs(b - c) / 2;
        bool condition3 = !mflag && std::abs(s - b) >= std::abs(c - d) / 2;
        bool condition4 = mflag && std::abs(b - c) < tolerance;
        bool condition5 = !mflag && std::abs(c - d) < tolerance;

        if (condition1 || condition2 || condition3 || condition4 || condition5) {
            // Use bisection
            s = (a + b) / 2;
            mflag = true;
        } else {
            mflag = false;
        }

        double fs = func(s);
        d = c;
        c = b;
        fc = fb;

        if (fa * fs < 0) {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }

        if (std::abs(fa) < std::abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }
    }

    return b;
}

// Root finding with multiple attempts for difficult cases.

// Correct solutions for Vm will range from ~ 10^-2 to 10^3
// (with the latter representing volumes for exterme T at low P)
// Often an incorrect solution can be found at Vm < 0.1 when the
// correct solution is in the order of 10^0, so try higher limits first
double findRoot(std::function<double(double)> func, double initialGuess,
                double tolerance = 1e-10) {

    // Try Brent's method with different bracketing intervals
    std::vector<std::pair<double, double>> brackets = {
                                                       {1.0, 1000.0},
                                                       {0.01, 10.0}, // be sure that we don't have a solution in this range before we look below 0.01
                                                       {0.001, 1.0}
                                                       //       {initialGuess * 0.1, initialGuess * 10},
                                                       //       {initialGuess * 0.01, initialGuess * 100},

                                                       };

    for (const auto& bracket : brackets) {
        double a = bracket.first;
        double b = bracket.second;

        try {
            double fa = func(a);
            double fb = func(b);

            // Check if we have a sign change (bracketing condition)
            if (fa * fb <= 0) {
                double result = brentMethod(func, a, b, tolerance);
                if (!std::isnan(result) && std::abs(func(result)) < 1e-6) {
                    return result;
                }
            }
        } catch (...) {
            continue; // Try next bracket
        }
    }

    return std::numeric_limits<double>::quiet_NaN();
}

// EOS function
double EOS_ZD09(double V, double Tm, double targetPm) {
    if (V <= 0) return std::numeric_limits<double>::infinity();

    double T2 = Tm * Tm;
    double T3 = T2 * Tm;
    double V2 = V * V;
    double V4 = V2 * V2;

    double Z = 1.0;
    Z += (c[0] + c[1]/T2 + c[2]/T3) / V;
    Z += (c[3] + c[4]/T2 + c[5]/T3) / V2;
    Z += (c[6] + c[7]/T2 + c[8]/T3) / V4;
    Z += (c[9] + c[10]/T2 + c[11]/T3) / (V4*V);
    Z += c[12]/T3/V2 * (c[13] + c[14]/V2) * std::exp(-c[14]/V2);

    return Rgas * Tm * Z / V - targetPm;
}

double ZD09volume(const std::array<double, 7>& vx, double Pbar, double TK, bool useK = true) {
    double epsilon = 0.0;
    double sigma = 0.0;

    for (int i = 0; i < 7; ++i) {
        for (int j = 0; j < 7; ++j) {
            double mk1 = useK ? K1[i][j] : 1.0;
            double mk2 = useK ? K2[i][j] : 1.0;
            epsilon += vx[i] * vx[j] * mk1 * std::sqrt(eps[i] * eps[j]);
            sigma += vx[i] * vx[j] * mk2 * (sig[i] + sig[j]) / 2.0;
        }
    }

    if (epsilon <= 0 || sigma <= 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // use epsilon(CH4)/sigma(CH4)^3 = 3.06258806161 rather than 3.0636
    double Pm = 3.06258806161 * std::pow(sigma, 3) * Pbar / epsilon;
    double Tm = 154.0 * TK / epsilon;

    // Define the function for root finding
    auto func = [Tm, Pm](double V) { return EOS_ZD09(V, Tm, Pm); };

    double Vm = findRoot(func, Vguess);

    if (std::isnan(Vm)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return 1000.0 * Vm * std::pow(sigma / 3.691, 3);
}

std::pair<double, std::array<double, 7>> ZD09fugacity(const std::array<double, 7>& vx,
                                                      double Pbar, double TK, bool useK = true) {
    double vol = ZD09volume(vx, Pbar, TK, useK) / 1000.0; // litre
    std::array<double, 7> fug = {};

    if (std::isnan(vol)) {
        return {0.0, fug};
    }

    double Z = Pbar * vol / Rgas / TK;
    std::array<double, 7> epsilons = {};
    std::array<double, 7> sigmas = {};
    double epsilon = 0.0;
    double sigma = 0.0;

    for (int i = 0; i < 7; ++i) {
        for (int j = 0; j < 7; ++j) {
            double mk1 = useK ? K1[i][j] : 1.0;
            double mk2 = useK ? K2[i][j] : 1.0;
            epsilons[i] += mk1 * vx[j] * std::sqrt(eps[i] * eps[j]);
            sigmas[i] += mk2 * vx[j] * (sig[i] + sig[j]) / 2.0;

            epsilon += vx[i] * vx[j] * mk1 * std::sqrt(eps[i] * eps[j]);
            sigma += vx[i] * vx[j] * mk2 * (sig[i] + sig[j]) / 2.0;
        }
    }

    //   double Pm = 3.0636 * std::pow(sigma, 3) * Pbar / epsilon;
    double Tm = 154.0 * TK / epsilon;
    double Vm = vol * std::pow(3.691 / sigma, 3);

    double S1 = (c[0] + c[1]/std::pow(Tm, 2) + c[2]/std::pow(Tm, 3)) / Vm;
    S1 += (c[3] + c[4]/std::pow(Tm, 2) + c[5]/std::pow(Tm, 3)) / (2.0 * std::pow(Vm, 2));
    S1 += (c[6] + c[7]/std::pow(Tm, 2) + c[8]/std::pow(Tm, 3)) / (4.0 * std::pow(Vm, 4));
    S1 += (c[9] + c[10]/std::pow(Tm, 2) + c[11]/std::pow(Tm, 3)) / (5.0 * std::pow(Vm, 5));
    S1 += c[12] / (2.0 * c[14] * std::pow(Tm, 3)) *
          (c[13] + 1.0 - (c[13] + 1.0 + c[14]/std::pow(Vm, 2)) * std::exp(-c[14]/std::pow(Vm, 2)));

    double S2 = (2.0*c[1]/std::pow(Tm, 2) + 3.0*c[2]/std::pow(Tm, 3)) / Vm;
    S2 += (2.0*c[4]/std::pow(Tm, 2) + 3.0*c[5]/std::pow(Tm, 3)) / (2.0 * std::pow(Vm, 2));
    S2 += (2.0*c[7]/std::pow(Tm, 2) + 3.0*c[8]/std::pow(Tm, 3)) / (4.0 * std::pow(Vm, 4));
    S2 += (2.0*c[10]/std::pow(Tm, 2) + 3.0*c[11]/std::pow(Tm, 3)) / (5.0 * std::pow(Vm, 5));
    S2 += 3.0 * c[12] / (2.0 * c[14] * std::pow(Tm, 3)) *
          (c[13] + 1.0 - (c[13] + 1.0 - c[14]/std::pow(Vm, 2)) * std::exp(-c[14]/std::pow(Vm, 2)));

    // Calculate fugacities
    for (int i = 0; i < 7; ++i) {
        double lnGamma = Z - 1.0 - std::log(Z) + S1 - 2.0*S2*(1.0 - epsilons[i]/epsilon) +
                         6.0 * (1.0 - Z) * (1.0 - sigmas[i]/sigma);
        fug[i] = Pbar * std::exp(lnGamma) * vx[i];
    }

    vol *= 1000.0;
    return {vol, fug};
}

class ZD09Window : public QMainWindow {
    Q_OBJECT

public:
    ZD09Window(QWidget *parent = nullptr) : QMainWindow(parent) {
        setupUI();
        connectSignals();
        setUnit(); // Initialize units
    }

private slots:
    void setUnit() {
        QString strP = pfactor->currentText();
        QString strT = toffset->currentText();
        setP->setSuffix("  " + strP);
        setP->setMaximum(2500000.0 / pfactor->currentData().toDouble());
        setT->setSuffix(" " + strT);
        flabel->setText("<SPAN STYLE=\"font-family:'Times New Roman'\"><i>f</i></SPAN><sub>i</sub> (" + strP + ")");
        calculate();
    }

    void clearText() {
        for (int i = 0; i < 7; ++i) {
            Fout[i]->setText("");
        }
        Vout->setText("");
    }

    void calculate() {
        std::array<double, 7> xx = {};

        // Get mole fractions
        for (int i = 0; i < 7; ++i) {
            xx[i] = Xin[i]->value();
        }

        // Normalize if sum != 1
        double sum = 0.0;
        for (int i = 0; i < 7; ++i) {
            sum += xx[i];
        }

        if (sum > 0 && std::abs(sum - 1.0) > 1e-10) {
            for (int i = 0; i < 7; ++i) {
                xx[i] /= sum;
                xx[i] = std::round(xx[i] * 100.0) / 100.0; // Round to 2 decimal places
                Xin[i]->blockSignals(true);
                Xin[i]->setValue(xx[i]);
                Xin[i]->blockSignals(false);
            }
        }

        if (sum == 0.0) return;

        double p = setP->value() * pfactor->currentData().toDouble();
        double t = setT->value() + toffset->currentData().toDouble();

        auto result = ZD09fugacity(xx, p, t);
        double volume = result.first;
        auto fugacities = result.second;

        // Convert fugacities to display units
        double pfac = pfactor->currentData().toDouble();
        for (int i = 0; i < 7; ++i) {
            fugacities[i] /= pfac;
            Fout[i]->setText(QString::number(fugacities[i], 'g', 7));
        }

        Vout->setText(QString::number(volume, 'g', 7));
    }

private:
    void setupUI() {
        setWindowTitle("Zhang and Duan 2009");

        auto centralWidget = new QWidget(this);
        setCentralWidget(centralWidget);

        auto mainLayout = new QVBoxLayout(centralWidget);
        auto layout = new QGridLayout();

        // Temperature offset combo
        toffset = new QComboBox();
        toffset->addItem("°C", 273.15);
        toffset->addItem("K", 0);

        // Pressure factor combo
        pfactor = new QComboBox();
        pfactor->addItem("Pa", 1.0e-5);
        pfactor->addItem("bar", 1.0);
        pfactor->addItem("MPa", 10.0);
        pfactor->addItem("kbar", 1000.0);
        pfactor->addItem("GPa", 1.0e+4);
        pfactor->setCurrentIndex(1);

        // Pressure spinbox
        setP = new QDoubleSpinBox();
        setP->setMinimum(0);
        setP->setMaximum(250);
        setP->setSuffix(" bar");
        setP->setValue(1);
        setP->setStepType(QAbstractSpinBox::AdaptiveDecimalStepType);
        setP->setAlignment(Qt::AlignCenter);

        // Temperature spinbox
        setT = new QDoubleSpinBox();
        setT->setMinimum(0);
        setT->setMaximum(5000);
        setT->setSuffix(" °C");
        setT->setValue(1000);
        setT->setStepType(QAbstractSpinBox::AdaptiveDecimalStepType);
        setT->setAlignment(Qt::AlignCenter);

        button = new QPushButton("Calculate");
        button->setMinimumSize(50,30);

        auto title = new QLabel("Calculate fluid volume and partial fugacities (<SPAN STYLE=\"font-family:'Times New Roman'\"><i>f</i></SPAN><sub>i</sub>) of fluid species using the equation of state of Zhang and Duan (2009)");
        title->setWordWrap(true);

        flabel = new QLabel("<SPAN STYLE=\"font-family:'Times New Roman'\"><i>f</i></SPAN><sub>i</sub> (bar)");
        flabel->setAlignment(Qt::AlignCenter);

        // Add widgets to layout
        layout->addWidget(title, 0, 0, 3, 4);
        layout->addWidget(pfactor, 3, 1);
        layout->addWidget(toffset, 4, 1);
        layout->addWidget(new QLabel("Pressure"), 3, 0, Qt::AlignRight);
        layout->addWidget(setP, 3, 2);
        layout->addWidget(new QLabel("Temperature"), 4, 0, Qt::AlignRight);
        layout->addWidget(setT, 4, 2);
        layout->addWidget(button, 5, 1, 1, 2);
        layout->addWidget(new QLabel("Species i"), 6, 0, Qt::AlignHCenter);
        layout->addWidget(new QLabel("X<sub>i</sub>"), 6, 1, Qt::AlignHCenter);
        layout->addWidget(flabel, 6, 2);

        // Create species input/output widgets
        for (int i = 0; i < 7; ++i) {
            Xin[i] = new QDoubleSpinBox();
            Xin[i]->setRange(0, 1);
            Xin[i]->setStepType(QAbstractSpinBox::AdaptiveDecimalStepType);
            Xin[i]->setAlignment(Qt::AlignCenter);

            Fout[i] = new QLineEdit("");
            Fout[i]->setAlignment(Qt::AlignRight);
            Fout[i]->setReadOnly(true);

            layout->addWidget(new QLabel(species[i]), 7 + i, 0, Qt::AlignHCenter);
            layout->addWidget(Xin[i], 7 + i, 1);
            layout->addWidget(Fout[i], 7 + i, 2);
        }

        Xin[1]->setValue(1.0); // Default to pure H2O

        // Separator line
        auto line = new QFrame();
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Plain);
        line->setLineWidth(1);
        layout->addWidget(line, 15, 0, 1, 3);

        // Volume output
        layout->addWidget(new QLabel("Volume (cc/mol):"), 16, 1);
        Vout = new QLineEdit("");
        Vout->setReadOnly(true);
        Vout->setAlignment(Qt::AlignRight);
        layout->addWidget(Vout, 16, 2);

        mainLayout->addLayout(layout);
    }

    void connectSignals() {
        connect(pfactor, &QComboBox::currentTextChanged, this, &ZD09Window::setUnit);
        connect(toffset, &QComboBox::currentTextChanged, this, &ZD09Window::setUnit);
        connect(button, &QPushButton::clicked, this, &ZD09Window::calculate);
        connect(setP, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &ZD09Window::calculate);
        connect(setT, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &ZD09Window::calculate);

        for (int i = 0; i < 7; ++i) {
            connect(Xin[i], QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &ZD09Window::clearText);
        }
    }

private:
    QComboBox *toffset, *pfactor;
    QDoubleSpinBox *setP, *setT;
    QPushButton *button;
    QLabel *flabel;
    std::array<QDoubleSpinBox*, 7> Xin;
    std::array<QLineEdit*, 7> Fout;
    QLineEdit *Vout;
};

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    // Initialize mixing parameters
    initializeMixingParameters();

    ZD09Window window;
    window.show();

    return app.exec();
}

#include "main.moc"
