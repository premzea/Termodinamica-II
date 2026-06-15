class Flash():

    def __init__(self, zi_F, temperature_f, pressure_f, TcDato_f, PcDato_f, wDato_f):
        self.zi = zi_F
        self.T = temperature_f
        self.P = pressure_f
        self.Tc = TcDato_f
        self.Pc = PcDato_f
        self.w = wDato_f

    def wilson(self):
        # Ecuación wilson
        lnKi = np.log(self.Pc / self.P) + 5.373 * (1 + self.w) * (1 - self.Tc / self.T)
        self.Ki = np.exp(lnKi)
        return self.Ki

    def beta(self):
        # Estimación de la fracción de fase de vapor en el sistema
        self.Ki = self.wilson()
        #Bmin = np.divide((self.Ki * self.zi - 1), (self.Ki - 1))
        Bmin = (self.Ki * self.zi - 1) / (self.Ki - 1)

        #print (("Bmin_inter = ", Bmin))

        Bmax = (1 - self.zi) / (1 - self.Ki)
        #print (("Bmax_inter = ", Bmax))
        self.Bini = (np.max(Bmin) + np.min(Bmax)) / 2
        print("inib =", self.Bini)
        return self.Bini

    def rice(self):
        # Ecuación de Rachford-Rice para el equilibrio líqudo-vapor
        self.fg = np.sum(self.zi * (self.Ki - 1) / (1 - self.Bini + self.Bini * self.Ki))
        self.dfg = - np.sum(self.zi * (self.Ki - 1) ** 2 / (1 - self.Bini + self.Bini * self.Ki) ** 2)
        #print g, dg
        return self.fg, self.dfg

    def composicion_xy(self):
        # Ecuación de Rachford-Rice para calcular la composición del equilibrio líqudo-vapor
        self.xi = self.zi / (1 - self.Bini + self.Bini * self.Ki)
        self.yi = (self.zi * self.Ki) / (1 - self.Bini + self.Bini * self.Ki)
        self.li = (self.zi * (1 - self.Bini)) / (1 - self.Bini + self.Bini * self.Ki)
        self.vi = (self.zi * self.Bini * self.Ki) / (1 - self.Bini + self.Bini * self.Ki)

        return self.xi, self.yi, self.li, self.vi

    def flash_ideal(self):
        # Solución del flash (T,P,ni) isotermico para Ki_(T,P)
        self.Bini = self.beta()
        self.Ki = self.wilson()
        # print ("Ki_(P, T) = ", self.Ki)
        Eg = self.rice()
        errorEq = abs(Eg[0])
        # Especificaciones del método Newton precario, mientras se cambia por una librería Scipy
        i, s, ep = 0, 1, 1e-5

        while errorEq > ep:
            Eg = self.rice()
            self.Bini = self.Bini - s * Eg[0] / Eg[1]
            errorEq = abs(Eg[0])
            i += 1
            if i >= 50:
                break

        xy = self.composicion_xy()
        print ("-"*53, "\n", "-"*18, "Mole fraction", "-"*18, "\n","-"*53)
        print ("\n", "-"*13, "Zi phase composition", "-"*13, "\n")
        print ("{0} = {1} \n {2} = {3} \n {4}={5} \n {6}={7} \n".format(Componentes_f1.value, self.zi[0], Componentes_f2.value, self.zi[1], Componentes_f3.value, self.zi[2], Componentes_f4.value, self.zi[3]))
        print ("Sumatoria zi = {0}".format(np.sum(self.zi)))
        print ("\n", "-"*13, "Liquid phase composition", "-"*13, "\n")
        print ("{0} = {1} \n {2} = {3} \n {4}={5} \n {6}={7} \n".format(Componentes_f1.value, self.xi[0], Componentes_f2.value, self.xi[1], Componentes_f3.value, self.xi[2], Componentes_f4.value, self.xi[3]))
        print ("Sumatoria xi = {0}".format(np.sum(self.xi)))
        print ("\n", "-"*14, "Vapor phase composition", "-"*13, "\n")
        print ("{0} = {1} \n {2} = {3} \n {4}={5} \n {6}={7} \n".format(Componentes_f1.value, self.yi[0], Componentes_f2.value, self.yi[1], Componentes_f3.value, self.yi[2], Componentes_f4.value, self.yi[3]))
        print ("Sumatoria yi = {0}".format(np.sum(self.yi)))
        print ("-"*53, "\n","Beta = {0}".format(self.Bini), "\n")
        print ("\n","Función R&R = {0}".format(Eg[0]), "\n")
        print ("\n","Derivada función R&R = {0}".format(Eg[1]), "\n", "-"*53)


        return #Eg[0], Eg[1], self.Bini

class FlashHP(Fugacidad, Flash):

    def __init__(self, zF):
            Fugacidad.__init__(self, eq, w, Tc, Pc, Tr, R, ep, ni, nT, nC, V, T, P, kij, lij, delta_1, k, Avsl)
            self.zF = zF



    def flash_PT(self):
            # Solución del flash (T,P,ni) isotermico para Ki_(T,P,ni)
            flashID = self.flash_ideal()
            print ("flash (P, T, zi)")
            print ("g, dg, B = ", flashID)
            print ("-"*66)

            self.Bini = flashID[2]
            print ("Beta_r ini = ", self.Bini)
            moles = self.composicion_xy()

            self.xi, self.yi = moles[0], moles[1]
            nil, niv = moles[2], moles[3]

            fi_F = self.fugac()

            self.Ki = fi_F[0] / fi_F[1]

            L = 1.0

            self.Ki = self.Ki * L

            Ki_1 = self.Ki
            print ("Ki_(P, T, ni) primera = ", self.Ki)

            print ("-"*66)

            #self.Ki = np.array([1.729, 0.832, 0.640])

            #self.Ki = self.wilson(self.Pc, self.Tc, self.w, self.T)
            #print "Ki_(P, T) = ", self.Ki

            while 1:
            i, s = 0, 0.1

            while 1:
                    Eg = self.rice()
                    print (Eg)
                    self.Bini = self.Bini - s * Eg[0] / Eg[1]
                    print (self.Bini)
                    errorEq = abs(Eg[0])
                    i += 1
                    #print i

                    #if self. Bini < 0 or self.Bini > 1:
                    #break
                    #    self.Bini = 0.5
                    if i >= 50:
                    pass
                    #break
                    if errorEq < 1e-5:
                    break

            print ("Resultado Real = ", Eg)
            print (" Beta r = ", self.Bini)

            moles = self.composicion_xy(zi, self.Ki, self.Bini)
            self.xi, self.yi = moles[0], moles[1]

            #xy = self.composicion_xy(zi, self.Ki, self.Bini)

            print ("C1 -i-C4 n-C4")
            print ("-"*13, "Composición de fase líquida", "-"*13)
            print ("xi = ", moles[0])
            print ("Sxi = ", np.sum(moles[0]))
            print ("-"*13, "Composición de fase vapor", "-"*13)
            print ("yi = ", moles[1])
            print ("Syi = ", np.sum(moles[1]))

            fi_F = self.fugac()

            self.Ki = fi_F[0] / fi_F[1]
            Ki_2 = self.Ki
            dKi = abs(Ki_1 - Ki_2)
            Ki_1 = Ki_2
            print ("Ki_(P, T, ni) = ", self.Ki)

            fun_Ki = np.sum(dKi)
            print ("fun_Ki = ", fun_Ki)

            if fun_Ki < 1e-5:
                    break

            return flashID