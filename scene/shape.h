//base copied from CIS 561

#pragma once


template<class T, int dim>
class Shape
{
public:
    Shape(std::string file) : filepath(file), isMoving(true) {}
    Shape() : filepath(""), isMoving(false) {}

    virtual ~Shape(){}
    virtual bool checkCollisions(const Eigen::Matrix<T, dim, 1> &pos, Eigen::Matrix<T, dim, 1> &out_pos) const = 0;
    void setCenter(Eigen::Matrix<T, dim, 1> &n_cen);
    void setVelocity(Eigen::Matrix<T, dim, 1> &n_vel);
    void outputFrame(int frame);
    void updatePosition(T dt);

protected:
    Eigen::Matrix<T, dim, 1> center; //Default initialize to zero
    Eigen::Matrix<T, dim, 1> velocity; //Default initialize to zero
    std::string filepath;
    bool isMoving;
    
};


template<class T, int dim>
void Shape<T, dim>::setCenter(Eigen::Matrix<T, dim, 1> &n_cen) {
    center = n_cen;
}

template<class T, int dim>
void Shape<T, dim>::setVelocity(Eigen::Matrix<T, dim, 1> &n_vel) {
    velocity = n_vel;
}

template<class T, int dim>
void Shape<T, dim>::updatePosition(T dt) {
    //center = center + (dt * velocity);
    Eigen::Matrix<T, dim, 1> temp = center + (dt * velocity);
    center = temp;
}

template<class T, int dim>
void Shape<T,dim>::outputFrame(int frame){
    if (isMoving) {
        Partio::ParticlesDataMutable* parts = Partio::create();
        Partio::ParticleAttribute posH, vH;
        posH = parts->addAttribute("position", Partio::VECTOR, 3);
        vH = parts->addAttribute("v", Partio::VECTOR, 3);

        int idx = parts->addParticle();
        float* p = parts->dataWrite<float>(posH, idx);
        float* v = parts->dataWrite<float>(vH, idx);
        for (int k = 0; k < 3; k++)
            p[k] = center[k];
        for (int k = 0; k < 3; k++)
            v[k] = velocity[k];

        // write frames to .bgeo file
        std::string f = std::to_string(frame);
        std::string particleFile = "";
        if(f.length() == 1)
            particleFile = "frame000" + f +".bgeo";
        else if(f.length() == 2)
            particleFile = "frame00" + f +".bgeo";
        else if(f.length() == 3)
            particleFile = "frame0" + f +".bgeo";
        else
            particleFile = "frame" + f +".bgeo";

        particleFile = filepath + "/" + particleFile;
        Partio::write(particleFile.c_str(), *parts);
        parts->release();
    }
}
